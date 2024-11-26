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

# Cargar el dataset
data = pd.read_csv("sequence_analysis_results.csv")

# Definir columnas no numéricas
non_numerical_columns = ['Type', 'Accession', 'UniProtID', 'SequenceLength']

# Eliminar columnas no numéricas para obtener datos numéricos
numerical_data = data.drop(columns=non_numerical_columns, errors='ignore')

# Manejar valores faltantes en los datos numéricos
numerical_data = numerical_data.fillna(numerical_data.mean())

# Estandarizar los datos
scaler = StandardScaler()
scaled_data = scaler.fit_transform(numerical_data)

# Convertir de nuevo a DataFrame para un manejo más fácil
scaled_df = pd.DataFrame(scaled_data, columns=numerical_data.columns)

# Aplicar k-Means clustering
kmeans = KMeans(n_clusters=2, random_state=42)
data['Cluster'] = kmeans.fit_predict(scaled_df)

# Reducir dimensiones para visualización con PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_df)

# Añadir resultados de PCA al DataFrame
data['PCA1'] = pca_result[:, 0]
data['PCA2'] = pca_result[:, 1]

# Visualizar los clusters
plt.figure(figsize=(8, 6))
sns.scatterplot(x='PCA1', y='PCA2', hue='Cluster', data=data, palette='viridis')
plt.title('Clusters Visualized with PCA')
plt.show()

# Información sobre los centroides
centroids = scaler.inverse_transform(kmeans.cluster_centers_)
centroids_df = pd.DataFrame(centroids, columns=numerical_data.columns)
print("Cluster Centroids:")
print(centroids_df)

# Dividir datos para entrenamiento y prueba
X = scaled_df
y = data['Cluster']

# Balancear las clases con SMOTE
smote = SMOTE(random_state=42)
X_resampled, y_resampled = smote.fit_resample(X, y)

# Dividir en entrenamiento y prueba
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=42)

# Sistema de épocas
max_epochs = 100
target_accuracy = 0.90  # Precisión objetivo
history = {"accuracy": [], "loss": []}  # Historial para graficar

# Configuración del modelo
mlp = MLPClassifier(
    hidden_layer_sizes=(100,),
    activation='relu',
    solver='adam',
    random_state=42,
    max_iter=1,  # Entrenaremos manualmente por épocas
    warm_start=True  # Permite continuar el entrenamiento sin reiniciar el modelo
)

# Entrenamiento por épocas
for epoch in range(1, max_epochs + 1):
    mlp.fit(X_train, y_train)  # Entrenamiento de una época

    # Predicciones en el conjunto de prueba
    y_pred = mlp.predict(X_test)

    # Evaluación del modelo
    epoch_accuracy = accuracy_score(y_test, y_pred)
    epoch_loss = 1 - epoch_accuracy  # Usamos pérdida como (1 - accuracy)

    # Almacenar historial
    history["accuracy"].append(epoch_accuracy)
    history["loss"].append(epoch_loss)

    print(f"Epoch {epoch}/{max_epochs} - Accuracy: {epoch_accuracy:.2f}, Loss: {epoch_loss:.2f}")

    # Detener el entrenamiento si alcanzamos la precisión objetivo
    if epoch_accuracy >= target_accuracy:
        print(f"\nTarget accuracy of {target_accuracy * 100}% reached at epoch {epoch}. Stopping training.")
        break

# Reporte final
print("\nFinal Classification Report:")
print(classification_report(y_test, y_pred))
print(f"Final Accuracy: {epoch_accuracy:.2f}")
print("Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))

# Gráficas de desempeño
plt.figure(figsize=(12, 5))

# Precisión por época
plt.subplot(1, 2, 1)
plt.plot(range(1, len(history["accuracy"]) + 1), history["accuracy"], label='Accuracy', color='blue')
plt.axhline(y=target_accuracy, color='red', linestyle='--', label='Target Accuracy')
plt.title('Accuracy per Epoch')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.legend()

# Pérdida por época
plt.subplot(1, 2, 2)
plt.plot(range(1, len(history["loss"]) + 1), history["loss"], label='Loss', color='orange')
plt.title('Loss per Epoch')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()

plt.tight_layout()

# Visualizar los centroides en el espacio 2D (usando PCA)
plt.figure(figsize=(8, 6))

# Reducir dimensiones para los centroides (para graficarlos en 2D)
pca_centroids = pca.transform(centroids)

# Plot de los clusters (agregar puntos de los datos)
sns.scatterplot(x='PCA1', y='PCA2', hue='Cluster', data=data, palette='viridis', alpha=0.6)

# Plot de los centroides
plt.scatter(pca_centroids[:, 0], pca_centroids[:, 1], s=200, c='red', marker='X', label='Centroids')

plt.title('Cluster Centroids')
plt.legend()
plt.show()

# Matriz de confusión
cm = confusion_matrix(y_test, y_pred)

# Crear el heatmap para visualizar la matriz de confusión
plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=['Cluster 0', 'Cluster 1'], yticklabels=['Cluster 0', 'Cluster 1'])
plt.title('Confusion Matrix')
plt.xlabel('Predicted')
plt.ylabel('True')
plt.show()
