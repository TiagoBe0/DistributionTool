import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Configuración de estilo
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_context("paper", font_scale=1.5)

# Nombres de las columnas
column_names = ['id', 'type', 'pc1', 'pc2', 'dist_to_ref', 'defect_type']

# Cargar datos: separador espacios, omitir primera línea (comentario)
df = pd.read_csv("pca_output.csv", sep=r'\s+', skiprows=1, names=column_names)

# Verificar que se cargó correctamente
print(df.head())
print(df.columns)

# Crear figura con dos subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Scatter plot PC1 vs PC2 coloreado por 'type'
scatter = ax1.scatter(df['pc1'], df['pc2'], c=df['type'], cmap='viridis',
                      s=20, alpha=0.7, edgecolors='none')
ax1.set_xlabel('PC1', fontsize=14)
ax1.set_ylabel('PC2', fontsize=14)
ax1.set_title('Componentes principales', fontsize=16)
ax1.grid(True, linestyle='--', alpha=0.5)
cbar = fig.colorbar(scatter, ax=ax1)
cbar.set_label('Tipo', fontsize=12)

# Histograma de dist_to_ref
sns.histplot(df['dist_to_ref'], bins=50, kde=True, ax=ax2, color='steelblue')
ax2.set_xlabel('Distancia a la referencia', fontsize=14)
ax2.set_ylabel('Frecuencia', fontsize=14)
ax2.set_title('Distribución de dist_to_ref', fontsize=16)
ax2.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('pca_plots.png', dpi=300, bbox_inches='tight')
plt.savefig('pca_plots.pdf', bbox_inches='tight')
plt.show()

# Boxplot por tipo
fig2, ax3 = plt.subplots(figsize=(10, 6))
sns.boxplot(data=df, x='type', y='dist_to_ref', palette='Set2', ax=ax3)
ax3.set_xlabel('Tipo', fontsize=14)
ax3.set_ylabel('Distancia a la referencia', fontsize=14)
ax3.set_title('Distancia por tipo', fontsize=16)
ax3.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('dist_by_type.png', dpi=300, bbox_inches='tight')
plt.savefig('dist_by_type.pdf', bbox_inches='tight')
plt.show()
