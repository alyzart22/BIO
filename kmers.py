#Programa para fragmentar en kmeros - Alida Zarate - 20 marzo 2021


import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO

class kmers():

	def __init__(self,long_kmer,saltos):
		self.long_kmer=long_kmer
		self.saltos=saltos


	#Obtiene del archivo fasta las secuencias y sus tamanos a un arreglo 
	def open_file(self):
		seq1 = SeqIO.parse(r"C:/Users/Alien/Documents/python/bio/sars_c.fasta", "fasta")
		registros=[]
		tamanos=[]
		for i in seq1:
  			ancho= len(i.seq)
  			tamanos.append(len(i.seq))
  			registros.append(i.seq)
		return registros, tamanos


	#Obtiene el total de kameros por secuencia
	def total_kmers_sec(self, x_tamano ):
		total_k=((x_tamano-self.long_kmer)+1) /self.saltos
		total_k=int(total_k)
		#print("---->", total_k)
		return total_k


	#Obtiene los kmeros en si y los guarda en una lista
	def proceso(self, kmeros, secuencia):
		sub_fragmentos=[]
		fragmentos=[]
		kmer_aux=[]
		fragmentos_aux=[]

		for i in range(kmeros):
			#print(i, "---")
			for k in range(self.long_kmer):
				sub_fragmentos.append(secuencia[k+(self.saltos*i)])

			kmer_aux=sub_fragmentos.copy()
			fragmentos.append(kmer_aux)
			sub_fragmentos.clear()
		#print(fragmentos)
		return fragmentos
		#print(fragmentos)




#----------------------------------------------------------------------------------------------------------

objeto=kmers(long_kmer = 20, saltos= 10)

secuencias, tamano_secuencias=objeto.open_file()
total=len(secuencias)
#print(total)
#Creamos la lista principal la cual tendra todos los kmers
all_fragmentos=[]

for i in range(1):
	kmeros=objeto.total_kmers_sec(tamano_secuencias[i])
	#print("#####################################################################")
	all_fragmentos.append(objeto.proceso(kmeros, secuencias[i]))
	#print(all_fragmentos)
print("---------------------------------------------")
print(all_fragmentos)