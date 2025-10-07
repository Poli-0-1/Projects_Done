#GENE-SET ENRICHMENT

library()
library()
library()
library()



#take results from other codes: p-values, log-fold change 
#Enrichment -> give a biological interpretetion to the data

#Many ways, i.e. : Gene Sets, Networks, Pathways

#Aim: break down cellular function in gene sets for insight gain
#Filter and understand which are the data of interest and that could be useful for the analysis, taking into account the statistical significance

#Most famous source for this is the GO(Gene Ontology)[GO POWER RANGERS!!!!!!!]

#GO -> frameworks used in bioinformatics, useful for finding a possible relationship between the data that we possess and the biological theory

#Many examples of different types of ontologies

#Presence of a hierarchy, imposing the GO as a DAG(Direct Acyclic Graph), going deeper in the hierarchy things  are more complex and specific

#Some GOs: 
#1)Molecular Function,describes activities taking place at molecular levels,
#not specifying the concept just the reaction and who takes part in it; 
#2)Biological Processes, try to capture a series of events resulting from 
#multiple ordered groups of molecular function, might be sometime difficult to distinguish the two of them
#but the typical difference is the presence of steps in the second one;
#3)Cellular Component, components in which a certain product can be found


#Transitive relationships in GO terms
#Part of relation captures "part of" the information, not certain information until we have data in hands

#Regulates, a process directly affects another one or its quality, this relationship one of the 
#elements may regulate another one, but at the same time in return there is a possibilit that it 
#will not be affected/regulated by the one that it affects

#GO needed in for avoiding language that may be ambigous, better orgnisation and 
#quicker/more efficient research 

#GO annotation: statement regarding a certain gene product with particular information
#included in the GEOterm

#Biological processes are not considered pathways, also due to the lack of wiring and "actors"

#MSigDB, db for groups of genes presented with similarities in a pathway with specific conditions/diseases

#It contains a collection of gene sets, could also be used for overlap checks with other products

#Methods for finding the enrichment:
          #Two-Class design: Expression matrix -> Ranked genes based on statistics -> Selective treshold 
          #Overlapping test between genes present in database and statistical-significant genes
          #Fisher test
          #Later focus on the score, depending on "Whole Distribution" or "Treshold-dependent"
          #Correlation to the phenotype
          #ES(Enrichment Score) checking level of the expression of the genes in the target, the higher the better
          #Estimate the empirical p-value by calculating all random gene sets the ES and search its position in the distribution that i calculated
          
