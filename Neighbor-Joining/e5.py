# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 16:42:46 2018

@author: alansoledar
"""
from Bio import Phylo
from io import StringIO

arq=open("Matrix.txt", "r")
lista = arq.readlines()
matrix = []
dic = {}
arvore = ""
matrix_len = len(lista)
for i in range(matrix_len):
    matrix.append([])

i=0
#preenche a matriz inicial
for x in lista:
    if x[-1] == '\n':
        x=x[:-1]
    matrix[i] = x.split(' ')
    i+=1
matrix[0].insert(0,'0')

while(matrix_len>3):
    #soma todos os S e coloca em um dicionário
    for y in range(1,matrix_len):
        soma = 0
        for x in range(1,matrix_len):
            soma += float(matrix[y][x])
        soma = soma/(matrix_len-3)
        dic[matrix[y][0]] = soma
        
    #encontra o menor valor M
    minimo = 1    
    for y in range(1,matrix_len):
        for x in range(1, y):
            soma = float(matrix[y][x])-dic[matrix[y][0]]-dic[matrix[x][0]]
        if soma < minimo:
            minimo = soma
            my = y
            mx = x
            
    #inicializaçaõ da lista U encabeçada pelo nome do ancestral comum que é a união dos 2 filhos e suas distancias do pai no formato [filho1:peso1, filho2:peso2]
    U =[str('('+matrix[my][0]+':'+str(round(((float(matrix[my][mx])/2)+(dic[matrix[my][0]]-dic[matrix[mx][0]])/2), 4))+', '+matrix[mx][0]+':'+str(round(((float(matrix[my][mx])/2)+(dic[matrix[mx][0]]-dic[matrix[my][0]])/2), 4))+')')]
    for y in range(1, matrix_len):
        U.append(str(round(((float(matrix[y][mx])+float(matrix[y][my])-float(matrix[mx][my]))/2), 4)))
            
    arvore = 2
    
    #adiciona uma nova linha usando a lista U como referência
    matrix = matrix+[U]
    #adiciona uma nova coluna usando a lista U como referência
    for i in range(len(U)):
        matrix[i].insert(matrix_len, U[i])
      
    #deleta as 2 linhas e 2 colunas referentes às especies q se juntaram
    del matrix[my]
    del matrix[mx]
    
    for y in range(matrix_len-1):
        del matrix[y][my]
        del matrix[y][mx]
    
    matrix_len-=1
    matrix[matrix_len-1].append('0')


#termina de montar a árvore juntando as ultimas 2 especies na matriz no formato de lista para printar com o biopython
treedata = '('+matrix[1][0]+':'+matrix[1][2]+', '+matrix[2][0]+')'
handle = StringIO(treedata)
tree = Phylo.read(handle, "newick")
#desenha a árvore
Phylo.draw(tree)
    
#printa a matriz final de forma rudimentar
for i in matrix:
    for j in i:
        print(j, end=' |\t')
    print('\n')
print(treedata)