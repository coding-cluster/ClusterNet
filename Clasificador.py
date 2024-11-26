import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix
from sklearn.decomposition import PCA
from imblearn.over_sampling import SMOTE
import matplotlib.pyplot as plt
import seaborn as sns

# Configuración general para visualización
sns.set(style="whitegrid")

# Cargar el dataset
data = pd.read_csv("sequence_analysis_results.csv")

# Definir columnas no numéricas
non_numerical_columns = ['Type', 'Accession', 'UniProtID', 'SequenceLength']

# Filtrar columnas no numéricas y trabajar solo con las numéricas
numerical_data = data.drop(columns=non_numerical_columns, errors='ignore')

# Manejar valores faltantes en los datos numéricos
numerical_data.fillna(numerical_data.mean(), inplace=True)

# Escalar los datos para normalización
scaler = StandardScaler()
scaled_data = scaler.fit_transform(numerical_data)

# Convertir de nuevo a DataFrame para manejo más fácil y consistencia
scaled_df = pd.DataFrame(scaled_data, columns=numerical_data.columns)

# Añadir columnas originales al DataFrame final
final_data = pd.concat([data[non_numerical_columns], scaled_df], axis=1)

# Aplicar k-Means clustering
kmeans = KMeans(n_clusters=3, random_state=42)
final_data['Cluster'] = kmeans.fit_predict(scaled_df)

# Reducir dimensiones para visualización con PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_df)

# Añadir resultados de PCA al DataFrame
final_data['PCA1'] = pca_result[:, 0]
final_data['PCA2'] = pca_result[:, 1]

# Visualización de clusters
plt.figure(figsize=(8, 6))
sns.scatterplot(x='PCA1', y='PCA2', hue='Cluster', data=final_data, palette='viridis')
plt.title('Clusters Visualized with PCA')
plt.show()

# Mostrar centroides del clustering
centroids = scaler.inverse_transform(kmeans.cluster_centers_)
centroids_df = pd.DataFrame(centroids, columns=numerical_data.columns)
print("Cluster Centroids:")
print(centroids_df)

# Preparar datos para clasificación
X = scaled_df
y = final_data['Cluster']

# Balancear las clases con SMOTE
smote = SMOTE(random_state=42)
X_resampled, y_resampled = smote.fit_resample(X, y)

# Dividir datos en entrenamiento y prueba
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=42)

# Configuración del modelo MLP
mlp = MLPClassifier(
    hidden_layer_sizes=(100,),
    activation='relu',
    solver='adam',
    random_state=42,
    max_iter=300  # Aumentado para asegurar buen entrenamiento
)

# Entrenamiento del modelo
mlp.fit(X_train, y_train)

# Evaluación del modelo
y_pred = mlp.predict(X_test)
print("\nFinal Classification Report:")
print(classification_report(y_test, y_pred))
print(f"Final Accuracy: {accuracy_score(y_test, y_pred):.2f}")
print("Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))

# Función de evaluación interactiva
def evaluate_model(features):
    """
    Evalúa el modelo dado un vector de características.
    """
    # Escalar las características
    features_scaled = scaler.transform([features])
    # Hacer la predicción
    prediction = mlp.predict(features_scaled)[0]
    # Mapear los clusters a tipos de virus
    cluster_map = {0: 'R5', 1: 'X4', 2: 'R5X4'}
    return cluster_map.get(prediction, "Unknown")

# Ejemplo de evaluación
example_features = scaled_df.iloc[0].tolist()  # Usar el primer registro como ejemplo
predicted_type = evaluate_model(example_features)
print(f"\nPredicted Type for Example Features: {predicted_type}")

# Matriz de confusión como heatmap
plt.figure(figsize=(8, 6))
cm = confusion_matrix(y_test, y_pred)
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=['R5', 'X4', 'R5X4'], yticklabels=['R5', 'X4', 'R5X4'])
plt.title('Confusion Matrix')
plt.xlabel('Predicted')
plt.ylabel('True')
plt.show()
