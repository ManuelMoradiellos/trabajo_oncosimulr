library(OncoSimulR)

#############################
###### PAPER PATHTIMEX ######
#############################

#### Ejemplo glioblastoma ####
## Fig. 3C, panel izquierdo

PT_glio <- allFitnessEffects(data.frame(
  parent = c("Root", "A", "Root", "Root", "B", "C", "D", "Root"),
  child = c("A", "B", "C", "D", "E", "E", "E", "F"),
  s = 0.2,
  sh = -1,
  typeDep = "OR"),
  geneToModule = c("Root" = "Root",
                   "A" = "CDKN2A, CDK4",
                   "B" = "TP53-M, MDM4, MDM2",
                   "C" = "NF1, FAF1, SPTA1, OBSCN, CNTNAP2",
                   "D" = "PTEN-M, PTEN-D, PIK3CA, IDH1",
                   "E" = "PDGFRA, LRP2",
                   "F" = "EGFR, RB1, TP53-D, PAOX"))
plot(PT_glio, expandModules = TRUE, autofit = TRUE)
## No funciona el evalAllGenotypes porque salen demasiados genotipos

#### Ejemplo cáncer colorrectal - datos Wood et al. 2007 - I ####
#### Esquema de la Fig. 3A, panel izquierdo

## Especificación del fitness y restricciones con módulos

PT_cr <- allFitnessEffects(data.frame(
  parent = c("Root", "Root", "Root", "A", "B", "C"),
  child = c("A", "C", "E", "B", "D", "D"),
  s = 0.5,
  sh = -1,
  typeDep = "MN"),
  geneToModule = c("Root" = "Root",
                   "A" = "APC",
                   "B" = "KRAS",
                   "C" = "TP53, EVC2",
                   "D" = "FBXW7, TCF7L2",
                   "E" = "PIK3CA, EPHA3"
  ))
plot(PT_cr, expandModules = TRUE)


#### Ejemplo cáncer colorrectal - datos Wood et al. 2007 - II ####
#### Siguiendo el esquema de la Fig. 4A, panel izquierdo

## Revisar diferencias con el ejemplo anterior
## Como son pocos genotipos, sí permite seguir adelante y hacer simulaciones

wood <- allFitnessEffects(data.frame(parent = c("Root", "A", "B", "C"),
                                     child = c("A", "B", "C", "D"),
                                     s = 0.1,
                                     sh = -0.01,
                                     typeDep = "MN"),
                          geneToModule = c("Root" = "Root",
                                           "A" = "APC, EPHA3",
                                           "B" = "EVC2, PIK3CA, TP53",
                                           "C" = "KRAS, TCF7L2",
                                           "D" = "FBXW7")
                          )

plot(wood, expandModules = TRUE, autofit = TRUE)

eag_wood <- evalAllGenotypes(wood)
plot(eag_wood, use_greppel = TRUE)

simuul_wood <- oncoSimulIndiv(wood,
                              model="McFL",
                              onlyCancer = FALSE,
                              finalTime = 10000,
                              verbosity = 0,
                              mu = 1e-5,
                              initSize = 5000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              detectionProb = NA,
                              detectionSize = NA,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
plot(simuul_wood, show = "genotypes", type = "line", 
     lwdClone = 2, legend.ncols = 4)

##############################
###### PAPER ¿OPCIONAL? ######
##############################

## Esquema de la Fig. 8C. Datos de glioblastoma multiforme
## Volvemos al problema de demasiados genotipos para trabajar:
## no podemos ir más allá de allFitnessEffects

njl <- allFitnessEffects(data.frame(parent = c("Root", "A", "PTEN", "TP53",
                                               "PIK3", "E", "F"),
                                    child = c("A", "PTEN", "TP53",
                                              "PIK3", "E", "F", "G"),
                                    s = 0.1,
                                    sh = -0.01,
                                    typeDep = "MN"),
                         geneToModule = c("Root" = "Root",
                                          "A" = "NF1, EGFR",
                                          "PTEN" = "PTEN",
                                          "TP53" = "TP53",
                                          "PIK3" = "PIK3R1, PIK3CA",
                                          "E" = "RB1, IDH1, STAG2",
                                          "F" = "ATRX, LZTR1",
                                          "G" = "BCOR, DCAF12L2"))

plot(njl, expandModules = TRUE, autofit = TRUE)

## Aquí ya se cuelga
eag_njl <- evalAllGenotypes(njl)
plot(eag_njl, use_greppel = TRUE)

## Si lo ejecutamos con tiempos bajos ok (no representativo), si no R aborta
simuul_njl <- oncoSimulIndiv(njl,
                             model="McFL",
                             onlyCancer = FALSE,
                             finalTime = 100,
                             verbosity = 0,
                             mu = 1e-5,
                             initSize = 5000,
                             keepPhylog = TRUE,
                             seed = NULL,
                             detectionProb = NA,
                             detectionSize = NA,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)
plot(simuul_njl, show = "genotypes", type = "line", lwdClone = 2)