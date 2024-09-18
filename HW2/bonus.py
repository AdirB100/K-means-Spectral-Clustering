import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import sklearn
from sklearn.cluster import KMeans
from sklearn import datasets
# allows the code to work on nova
matplotlib.use('Agg')

iris_dataset = datasets.load_iris(return_X_y=True)
iris_dataset = np.array(iris_dataset[0])
iris_dataset = pd.DataFrame(iris_dataset)
inertias = []

for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, init='k-means++', n_init=1, random_state=0)
    kmeans.fit(iris_dataset)
    inertias.append(kmeans.inertia_)

plt.plot(range(1, 11), inertias)
plt.title('Elbow Method for Selection of Optimal \"K\" clusters')
plt.xlabel('k')
plt.ylabel('Average Dispersion')
plt.annotate("Elbow Point", (3, inertias[2]), (5, 400), arrowprops=dict(facecolor='black', shrink=0.05),color='black')
plt.savefig("elbow.png")

