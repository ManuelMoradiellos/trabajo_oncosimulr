library(OncoSimulR)

##############################
#### PAPER CIRIELLO ET AL.####
##############################

## Empezamos aquí: no podemos hacer esto
m1 <- allFitnessEffects(data.frame(
  parent = c("Root", "Root"),
  child = c("A", "B"),
  s = 0.2,
  sh = -1,
  typeDep = "OR"),
  geneToModule = c("Root" = "Root",
                   "A" = "CDKN2A, CDK4, RB1",
                   "B" = "CDKN2B, CDK4, RB1")
  )
## ¿Es esto lo que realmente queremos hacer? (Problema imposible de resolver...)

#############################
#### POSIBLES SOLUCIONES ####
#############################

## Ejemplo Víctor especificando los fitnesses manualmente

## Problema: la transición CDK4/RB1 a CDKN2A/B + CD4/KRB1 
## disminuye el fitness --> ¿sentido biológico?
## Aparecen genotipos "prohibidos" en partes altas del fitness landscape

## Ventaja: parece que en las simulaciones ganan genotipos mutually exclusive

manual <- data.frame(Genotype = c("WT", "CDKN2A", "CDKN2B", "CDK4", "RB1", 
                                  "CDKN2A, CDKN2B", "CDKN2A, CDK4", "CDKN2A, RB1",
                                  "CDKN2B, CDK4", "CDKN2B, RB1", "RB1, CDK4",
                                  "CDKN2A, CDKN2B, CDK4", "CDKN2A, CDKN2B, RB1",
                                  "CDKN2A, CDK4, RB1", "CDKN2B, CDK4, RB1",
                                  "CDKN2A, CDKN2B, CDK4, RB1"), 
                     Fitness = c(1, 1.25, 1.25, 1.5, 1.5,
                                 1.5, 1.25, 1.25,
                                 1.25, 1.25, 1.5,
                                 1.5, 1.5,
                                 1.25, 1.25,
                                 1.5)
                    ) # Duda fitness CDKN2A + CDKN2B = 1.5 o 1.56 (= 1.25 * 1.25??

manual_aFE <- allFitnessEffects(genotFitness = manual)
plot(evalAllGenotypes(manual_aFE, order = FALSE, addwt = TRUE, ),
     use_ggrepel = TRUE)

simuul <- oncoSimulIndiv(manual_aFE,
                         model="McFL",
                         onlyCancer = FALSE,
                         finalTime = 1000,
                         verbosity = 0,
                         mu = 1e-5,
                         initSize = 5000,
                         keepPhylog = TRUE,
                         seed = NULL,
                         detectionProb = NA,
                         detectionSize = NA,
                         errorHitMaxTries = FALSE,
                         errorHitWallTime = FALSE)
plot(simuul,  show="genotypes", type="line", lwdClone = 2)


## Una alternativa: todos los fitnesses iguales, efecto multiplicativo,
## pero si tenemos CDK4/RB1 + CDKN2A/B y se añade el otro CDKN2 sube el fitness
## evitamos problema de que baje, pero incumple mutual exclusivity

## Sigue habiendo genotipos "prohibidos" en la parte alta del landscape

## Simulaciones: parece que cumplen mutual exclusivity (revisar)

manual_b <- data.frame(Genotype = c("WT", "CDKN2A", "CDKN2B", "CDK4", "RB1", 
                                  "CDKN2A, CDKN2B", "CDKN2A, CDK4", "CDKN2A, RB1",
                                  "CDKN2B, CDK4", "CDKN2B, RB1", "RB1, CDK4",
                                  "CDKN2A, CDKN2B, CDK4", "CDKN2A, CDKN2B, RB1",
                                  "CDKN2A, CDK4, RB1", "CDKN2B, CDK4, RB1",
                                  "CDKN2A, CDKN2B, CDK4, RB1"), 
                     Fitness = c(1, 1.25, 1.25, 1.5, 1.5,
                                 1.25*1.25, 1.25, 1.25,
                                 1.25, 1.25, 1.25,
                                 1.25*1.25, 1.25*1.25,
                                 1.25, 1.25,
                                 1.25*1.25)
                     )
manual_aFE_b <- allFitnessEffects(genotFitness = manual_b)
plot(manual_aFE_b)

eAG_b <- evalAllGenotypes(manual_aFE_b, order = FALSE, addwt = TRUE)
plot(eAG_b)

simuul_b <- oncoSimulIndiv(manual_aFE_b,
                         model="McFL",
                         onlyCancer = FALSE,
                         finalTime = 1000,
                         verbosity = 0,
                         mu = 1e-5,
                         initSize = 5000,
                         keepPhylog = TRUE,
                         seed = NULL,
                         detectionProb = NA,
                         detectionSize = NA,
                         errorHitMaxTries = FALSE,
                         errorHitWallTime = FALSE)
plot(simuul_b,  show="genotypes", type="line", lwdClone = 2)


## Última propuesta: poner una penalización cuando no se cumplen condiciones 
## de mutual exclusivity. En este ejemplo se multiplica por 0.1

## El fitness landscape evita genotipos prohibidos y parece que simulaciones tb.

## Problema: creemos que es una solución bastante artificial, el fitness baja 
## mucho al añadir mutaciones (no es letalidad sintética, en realidad "no sube",
## creemos que no tiene por qué bajar por debajo del wt)

## Simulaciones evitan mutual exclusivity (lógico porque el fitness baja mucho
## cuando nos desviamos)

manual_penalty <- data.frame(Genotype = c("WT", "CDKN2A", "CDKN2B", "CDK4", "RB1", 
                                      "CDKN2A, CDKN2B", "CDKN2A, CDK4", "CDKN2A, RB1",
                                      "CDKN2B, CDK4", "CDKN2B, RB1", "RB1, CDK4",
                                      "CDKN2A, CDKN2B, CDK4", "CDKN2A, CDKN2B, RB1",
                                      "CDKN2A, CDK4, RB1", "CDKN2B, CDK4, RB1",
                                      "CDKN2A, CDKN2B, CDK4, RB1"), 
                         Fitness = c(1, 1.25, 1.25, 1.56, 1.56,
                                     1.25*1.25, 1.25*0.1, 1.25*0.1,
                                     1.25*0.1, 1.25*0.1, 1.25*0.1,
                                     1.25*1.25*0.1, 1.25*1.25*0.1,
                                     1.25*0.1*0.1, 1.25*0.1*0.1,
                                     1.25*1.25*0.1*0.1)
)

manual_aFE_penalty <- allFitnessEffects(genotFitness = manual_penalty)
plot(evalAllGenotypes(manual_aFE_penalty, order = FALSE, addwt = TRUE, ), use_ggrepel = TRUE)

simuul_penalty <- oncoSimulIndiv(manual_aFE_penalty,
                         model="McFL",
                         onlyCancer = FALSE,
                         finalTime = 1000,
                         verbosity = 0,
                         mu = 1e-5,
                         initSize = 5000,
                         keepPhylog = TRUE,
                         seed = NULL,
                         detectionProb = NA,
                         detectionSize = NA,
                         errorHitMaxTries = FALSE,
                         errorHitWallTime = FALSE)
plot(simuul_penalty,  show="genotypes", type="line", lwdClone = 2)