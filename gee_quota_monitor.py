"""
Google Earth Engine - Monthly EECU Quota Monitor
=================================================
Checks your project's EECU-hour consumption for the current month
and compares it against the noncommercial tier limits.

Noncommercial Tier Limits (as of April 27, 2026):
  - Community Tier  :    150 EECU-hours / month
  - Contributor Tier:  1,000 EECU-hours / month
  - Partner Tier    : 100,000 EECU-hours / month

Requirements:
  pip install google-cloud-monitoring google-auth

Authentication:
  Run:  gcloud auth application-default login
  Then: gcloud auth application-default set-quota-project YOUR_PROJECT_ID
  (This silences the "no quota project" UserWarning)

Usage:
  python gee_quota_monitor.py --project YOUR_PROJECT_ID [--tier contributor]
"""

import argparse
import datetime
import sys
import warnings

# Suppress the "no quota project" ADC warning — we set it explicitly below
warnings.filterwarnings("ignore", message=".*quota project.*", category=UserWarning)

# ── Tier definitions ──────────────────────────────────────────────────────────
TIER_LIMITS_EECU_HOURS = {
    "community":   150,
    "contributor": 1_000,
    "partner":     100_000,
}

# ── The correct Cloud Monitoring metric (confirmed from Google's official notebook)
# earthengine.googleapis.com/project/cpu/usage_time  (EECU-seconds, label: compute_type)
METRIC_TYPE = "earthengine.googleapis.com/project/cpu/usage_time"

# ── Helpers ───────────────────────────────────────────────────────────────────

def get_month_window():
    """Return (start, end) UTC datetimes for the current calendar month."""
    now = datetime.datetime.now(datetime.timezone.utc)
    start = now.replace(day=1, hour=0, minute=0, second=0, microsecond=0)
    return start, now


def build_credentials(project_id: str):
    """
    Build ADC credentials and attach the given project as quota project.
    This avoids the 'no quota project' UserWarning.
    """
    try:
        import google.auth
        import google.auth.transport.requests
        from google.auth import impersonated_credentials  # noqa – just ensure import works
    except ImportError:
        sys.exit("Missing dependency: pip install google-auth")

    creds, detected_project = google.auth.default(
        scopes=["https://www.googleapis.com/auth/cloud-platform"]
    )

    # Attach quota project so Cloud Monitoring billing works correctly
    if hasattr(creds, "with_quota_project"):
        creds = creds.with_quota_project(project_id)

    return creds


def list_available_ee_metrics(project_id: str, creds) -> list:
    """
    Diagnostic helper: list all EE metric types visible for this project.
    Useful when the main metric returns 404.
    """
    try:
        from google.cloud import monitoring_v3
    except ImportError:
        return []

    client = monitoring_v3.MetricServiceClient(credentials=creds)
    resource_name = f"projects/{project_id}"

    try:
        descriptors = client.list_metric_descriptors(
            name=resource_name,
            filter='metric.type = starts_with("earthengine.googleapis.com")',
        )
        return [d.type for d in descriptors]
    except Exception:
        return []


def fetch_eecu_seconds(project_id: str, start: datetime.datetime,
                       end: datetime.datetime, creds) -> dict:
    """
    Query Cloud Monitoring for EE CPU usage and return EECU-seconds
    broken down by compute_type (online / batch).
    """
    try:
        from google.cloud import monitoring_v3
    except ImportError:
        sys.exit("Missing dependency: pip install google-cloud-monitoring")

    client = monitoring_v3.MetricServiceClient(credentials=creds)
    resource_name = f"projects/{project_id}"

    from google.protobuf.timestamp_pb2 import Timestamp
    start_ts = Timestamp(); start_ts.FromDatetime(start)
    end_ts   = Timestamp(); end_ts.FromDatetime(end)

    interval = monitoring_v3.TimeInterval(start_time=start_ts, end_time=end_ts)

    window_seconds = int((end - start).total_seconds())

    aggregation = monitoring_v3.Aggregation(
        alignment_period={"seconds": window_seconds},
        per_series_aligner=monitoring_v3.Aggregation.Aligner.ALIGN_SUM,
        cross_series_reducer=monitoring_v3.Aggregation.Reducer.REDUCE_SUM,
        group_by_fields=["metric.labels.compute_type"],
    )

    request = monitoring_v3.ListTimeSeriesRequest(
        name=resource_name,
        filter=f'metric.type = "{METRIC_TYPE}"',
        interval=interval,
        aggregation=aggregation,
        view=monitoring_v3.ListTimeSeriesRequest.TimeSeriesView.FULL,
    )

    results = {"online_eecu_seconds": 0.0, "batch_eecu_seconds": 0.0}

    try:
        for series in client.list_time_series(request=request):
            compute_type = series.metric.labels.get("compute_type", "unknown").lower()
            series_total = sum(p.value.double_value for p in series.points)
            if "online" in compute_type:
                results["online_eecu_seconds"] += series_total
            elif "batch" in compute_type:
                results["batch_eecu_seconds"] += series_total
            else:
                key = f"{compute_type}_eecu_seconds"
                results[key] = results.get(key, 0.0) + series_total

    except Exception as exc:
        err_str = str(exc)

        # ── 404: metric not found ────────────────────────────────────────────
        if "404" in err_str:
            print(f"\n⚠️  Metric '{METRIC_TYPE}' not found for project '{project_id}'.")
            print("   Possible reasons:")
            print("   1. Zero Earth Engine compute has run this month — the metric only")
            print("      appears after your first EE job completes in a billing period.")
            print("   2. The Cloud Monitoring API is not enabled. Enable it at:")
            print(f"      https://console.cloud.google.com/apis/library/monitoring.googleapis.com?project={project_id}")
            print("   3. Your account lacks the 'monitoring.viewer' IAM role on this project.")
            print("\n   Checking which EE metrics ARE available for your project …\n")
            available = list_available_ee_metrics(project_id, creds)
            if available:
                print("   Available EE metrics:")
                for m in available:
                    print(f"     • {m}")
            else:
                print("   No EE metrics found at all — the project may have no EE activity")
                print("   or the Cloud Monitoring API may be disabled.")
            print()
            sys.exit(1)

        # ── 403: permission denied ───────────────────────────────────────────
        elif "403" in err_str or "PERMISSION_DENIED" in err_str:
            sys.exit(
                f"\n❌ Permission denied: {exc}\n\n"
                f"   Grant yourself the 'Monitoring Viewer' role:\n"
                f"   https://console.cloud.google.com/iam-admin/iam?project={project_id}\n"
            )

        # ── Other errors ─────────────────────────────────────────────────────
        else:
            sys.exit(f"\n❌ Cloud Monitoring error: {exc}\n")

    results["total_eecu_seconds"] = (
        results["online_eecu_seconds"] + results["batch_eecu_seconds"]
    )
    return results


def s_to_h(s: float) -> float:
    return s / 3600.0


def format_bar(used: float, limit: float, width: int = 40) -> str:
    pct   = min(used / limit, 1.0) if limit > 0 else 0
    filled = int(pct * width)
    return f"[{'█' * filled}{'░' * (width - filled)}] {pct * 100:.1f}%"


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Monitor Google Earth Engine monthly EECU quota consumption."
    )
    parser.add_argument("--project", "-p", required=True,
                        help="Google Cloud project ID linked to your EE project.")
    parser.add_argument("--tier", "-t", default="contributor",
                        choices=list(TIER_LIMITS_EECU_HOURS.keys()),
                        help="Your noncommercial tier (default: contributor).")
    args = parser.parse_args()

    project_id      = args.project
    tier            = args.tier.lower()
    tier_limit_hours = TIER_LIMITS_EECU_HOURS[tier]

    start, end = get_month_window()

    print("=" * 62)
    print("  Google Earth Engine — Monthly EECU Quota Monitor")
    print("=" * 62)
    print(f"  Project  : {project_id}")
    print(f"  Tier     : {tier.capitalize()} ({tier_limit_hours:,} EECU-hours/month)")
    print(f"  Period   : {start.strftime('%Y-%m-%d')} → {end.strftime('%Y-%m-%d %H:%M')} UTC")
    print("-" * 62)

    print("\n⏳ Authenticating …")
    creds = build_credentials(project_id)

    print("⏳ Fetching usage data from Cloud Monitoring …\n")
    data = fetch_eecu_seconds(project_id, start, end, creds)

    online_h = s_to_h(data["online_eecu_seconds"])
    batch_h  = s_to_h(data["batch_eecu_seconds"])
    total_h  = s_to_h(data["total_eecu_seconds"])

    # Print any extra compute types (e.g. "interactive")
    extra_keys = [k for k in data if k not in
                  ("online_eecu_seconds", "batch_eecu_seconds", "total_eecu_seconds")]

    print(f"  Online (Code Editor / Python API)  : {online_h:>10.3f} EECU-hours")
    print(f"  Batch  (Export tasks)              : {batch_h:>10.3f} EECU-hours")
    for k in extra_keys:
        label = k.replace("_eecu_seconds", "").capitalize()
        print(f"  {label:<35}: {s_to_h(data[k]):>10.3f} EECU-hours")
    print(f"  {'─' * 53}")
    print(f"  Total consumed this month          : {total_h:>10.3f} EECU-hours")
    print(f"  Tier limit ({tier.capitalize():<12})          : {tier_limit_hours:>10,} EECU-hours")
    remaining = tier_limit_hours - total_h
    print(f"  Remaining                          : {max(remaining, 0):>10.3f} EECU-hours")

    print()
    print(f"  Usage : {format_bar(total_h, tier_limit_hours)}")
    print()

    pct = (total_h / tier_limit_hours * 100) if tier_limit_hours else 0

    if total_h > tier_limit_hours:
        over = total_h - tier_limit_hours
        print("🔴  STATUS: TIER EXCEEDED")
        print(f"    You have used {total_h:.2f} EECU-hours — {over:.2f} hours OVER the limit.")
        print("    Your project is likely in Restricted mode (slowed computations).")
        print("    → Consider upgrading to Partner tier or optimising your workflows.")
    elif pct >= 90:
        print("🟠  STATUS: CRITICAL — approaching limit")
        print(f"    {pct:.1f}% used. Only {remaining:.2f} EECU-hours remain.")
        print("    → Monitor closely to avoid hitting Restricted mode.")
    elif pct >= 70:
        print("🟡  STATUS: WARNING — moderate usage")
        print(f"    {pct:.1f}% used. {remaining:.2f} EECU-hours remaining.")
    else:
        print("🟢  STATUS: OK — within tier limit")
        print(f"    {pct:.1f}% used. {remaining:.2f} EECU-hours remaining.")

    print()
    print("-" * 62)
    print("  Tip: use ee.profilePrinting() (Python) or 'Run with Profiler'")
    print("  (Code Editor) to find expensive operations before they run.")
    print("  Docs: https://developers.google.com/earth-engine/guides/monitoring_usage")
    print("=" * 62)


if __name__ == "__main__":
    main()
