# -*- coding: utf-8 -*-
"""
Created on Sat Aug 28 18:20:44 2021

@author: diegogotex

"""
#bibliotecas
import argparse
from Bio import SeqIO
import pandas as pd
from pyfaidx import Fasta

###########
##PARTE 1##
###########

#definindo os arquivos
#construindo o passador de argumento
ap  = argparse.ArgumentParser()

#acicionando os argumentos
ap.add_argument('-f','-fasta', required=True, help="arquivo fasta")

ap.add_argument('-d','-depth_file', required=True, help="output from samtools depth")

args = vars(ap.parse_args())

#caminho para o fasta = args['f']
#caminho para o samtools depth = args['d']


###########
##PARTE 2##
###########

#pegar o cabecalho do fasta
for seq in SeqIO.parse(args['f'], "fasta"):
   header = seq.id

'''
o arquivo de profundidade de saida do samtools e assim:

samtools depth teste.bam > teste.DEPTH.txt    

Chr(REF) posicao cobertura
teste 1 5
teste 2 12
teste 3 13
teste 4 13
teste 5 6
teste 6 30
teste 7 31
teste 8 31
teste 9 3
teste 10 32
teste 11 11
teste 12 13
teste 13 13
teste 14 14
teste 15 12
teste 16 13
'''

###########
##PARTE 3##
###########

#olha o arquivo de cobertura e identifica quais bases estao abaixo do ponto de corte
#guarda a linha e adiciona um 'N' na ultima coluna
dp_file = pd.read_csv(args['d'], delimiter='\t', names=['chr', 'pos','depth'])

#criando uma tabela pra ser populada
low_cov_table = pd.DataFrame(columns=['header','pos','base'])

#loop para popular a tabela com as posicoes de baixa cobertura
for i in range(len(dp_file)):
    if dp_file.loc[i,"depth"] < 1:
        #junta o cabecalho da sequencia, posicao com cobertura abaixo de 10 e "N" em uma nova linha do df
        low_cov_table = low_cov_table.append({'header': header, 'pos': dp_file.loc[i,"pos"], 'base':'N'}, ignore_index=True)
        #print(header, dp_file.loc[i,"pos"], 'N', sep='\t')

#escrevendo o arquivo tabulado com as posicoes
low_cov_table.to_csv('low_cov_table.txt', sep='\t', header=False, index=False)


'''
o fasta deve ser assim:
    
>teste
ACTGACTGACTGACTG 


o arquivo de mutacao deve ser assim:     
    
teste 1 N
teste 5 N
teste 9 N
    
'''

###########
##PARTE 4##
###########


#abrindo o arquivo com as alteracoes
with open('low_cov_table.txt') as mut_table:
  # mutable Fasta modifica o arquivo de entrada
  # por tanto tenha certeza de que esta colocando o arquivo certo
  with Fasta(args['f'], mutable=True) as fasta:
    #lendo linha por linha do arquivo com as alteracoes
     for line in mut_table:
      #definindo o que significa cada coluna do arquivo
      rname, pos, base = line.rstrip().split()
      # convertendo as coordenadas baseadas em 1 para 0 
      fasta[rname][int(pos) - 1] = base
      
      
'''
a saida e assim: 

>teste
NCTGNCTGNCTGACTG   

REFs
https://www.datacamp.com/community/tutorials/argument-parsing-in-python
https://www.biostars.org/p/265811/ (matt shirley)
'''