import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Función para calcular el contenido GC de una secuencia
def GC_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)

## listas
GC_contents = []
GC_aleatorio = []
gene_positions = []

## Archivo FASTA y secuencias
with open("mycoplasma_genes.txt") as file:
    title_gene = None   #encabezado de cada gen
    current_sequence = ""
    counter = 0

    for line in file:
        if line.startswith(">"):
            if title_gene:
                GC = GC_content(current_sequence)
                GC_contents.append(GC)
                counter += 1
                gene_positions.append(counter)
                GC_aleatorio.append(GC)

            title_gene = line.strip()  # Guarda el encabezado del gen
            current_sequence = ""
        else:
            current_sequence += line.strip()  # Concatena la secuencia


# Reordenar aleatoriamente el contenido GC
random.shuffle(GC_aleatorio)


#  listas
#print(f"Contenido GC: {GC_contents}")

#print(f"\nPosición de cada gen: {gene_positions}")

#print(f"Contenido GC aleatorio: {GC_aleatorio}")

#Generar la tabla

data = {"Gene Position": gene_positions,
    "GC Content": GC_contents,
    "GC Content Aleatorio": GC_aleatorio}

df = pd.DataFrame(data)

# Mostrar el DataFrame en una tabla
#print(df)

#Valores anormales

umbral_superior = 0.40
puntos_anormales = df[df["GC Content"] > umbral_superior]
#print(puntos_anormales)


#Crear un lineplot
plt.figure(figsize=(16, 10))
sns.lineplot(data=df, x="Gene Position", y="GC Content", label="GC Content")
sns.lineplot(data=df, x="Gene Position", y="GC Content Aleatorio", label="GC Content Aleatorio")
plt.xlabel("Gene Position")
plt.ylabel("GC Content")
plt.title("Contenido GC vs. Gene Position")

for index, row in puntos_anormales.iterrows():     #Datos de esa fila
    plt.annotate(row["Gene Position"],
                 (row["Gene Position"],     #Posicion de etiqueta
                  row["GC Content"]),
                 textcoords="offset points",
                 xytext=(0, 10),
                 ha='center'
                 )
plt.savefig("grafica_mycoplasma.jpg", dpi=600)
plt.show()





