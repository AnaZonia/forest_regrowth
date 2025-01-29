import numpy as np

# Dados fornecidos
saldo_devedor = 55405
juros_anual = 0.05
juros_mensal = juros_anual / 12
amortizacao_mensal = 196
seguro_mensal = 30
correcao_mensal = 39
amortizacao_adicional = 0
investimento_dolar = 500
crescimento_dolar_mensal = 1+(0.0855/12)

# Variáveis para amortização do financiamento
meses = 0
juros_totais = 0
total_guardado_dolar = 0

# Loop para amortização
# while saldo_devedor > 0:
while meses < 24:
    meses += 1

    juros_mes = saldo_devedor * juros_mensal
    juros_totais += juros_mes

    # Atualizar retorno do investimento em dólar
    total_guardado_dolar += investimento_dolar 
    total_guardado_dolar *= crescimento_dolar_mensal

    # Atualizar saldo devedor
    saldo_devedor -= (amortizacao_mensal + amortizacao_adicional)

anos = meses // 12
meses_restantes = meses % 12

# Comparação
print(f"Tempo para quitar o financiamento: {anos} anos e {meses_restantes} meses")
print(f"Total de juros pagos: R${juros_totais:.2f}")
print(f"Total guardado no dólar: R${total_guardado_dolar:.2f}")
print(f"Retorno estimado no dólar: R${total_guardado_dolar-(meses*investimento_dolar):.2f}")
