import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import numpy as np

# 1. Cargar el archivo sin encabezados y asignar nombres de columnas
column_names = ['id', 'type', 'x', 'y', 'z', 'dist_to_ref', 'defect_prob', 'defect_type']
df = pd.read_csv('output.csv', sep=r'\s+', comment='#', header=None, names=column_names)

print(df.head())
print(f"Total de registros: {len(df)}")

# 2. Muestrear los datos para mejorar rendimiento (por ejemplo, 20000 puntos)
sample_size = 20000
if len(df) > sample_size:
    df_sample = df.sample(n=sample_size, random_state=42)
else:
    df_sample = df
print(f"Usando {len(df_sample)} puntos para los gráficos.")

# 3. Gráfico 3D interactivo con plotly (muestreo)
fig = px.scatter_3d(df_sample, x='x', y='y', z='z', color='defect_prob',
                     title='Defectos 3D - Probabilidad (muestreo)',
                     labels={'defect_prob': 'Probabilidad'},
                     opacity=0.6, size_max=3)
fig.show()

# 4. Gráfico 3D por tipo de defecto
fig2 = px.scatter_3d(df_sample, x='x', y='y', z='z', color='defect_type',
                     title='Defectos 3D - Tipo de defecto (muestreo)',
                     category_orders={'defect_type': sorted(df_sample['defect_type'].unique())})
fig2.show()

# 5. Proyección XY con matplotlib (puedes usar todos los datos si quieres)
plt.figure(figsize=(10,8))
sc = plt.scatter(df['x'], df['y'], c=df['defect_prob'], cmap='viridis', s=1, alpha=0.5)
plt.colorbar(sc, label='Probabilidad de defecto')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Proyección XY de defectos (todos los datos)')
plt.show()

# 6. Histograma de probabilidad
plt.figure(figsize=(10,6))
plt.hist(df['defect_prob'], bins=100, edgecolor='black')
plt.xlabel('Probabilidad de defecto')
plt.ylabel('Frecuencia')
plt.title('Distribución de defect_prob')
plt.show()

# 7. Conteo por tipo de defecto
defect_counts = df['defect_type'].value_counts()
plt.figure(figsize=(10,6))
defect_counts.plot(kind='bar')
plt.xlabel('Tipo de defecto')
plt.ylabel('Cantidad')
plt.title('Cantidad de defectos por tipo')
plt.xticks(rotation=45)
plt.show()
