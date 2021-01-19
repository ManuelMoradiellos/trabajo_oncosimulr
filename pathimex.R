library(OncoSimulR)
set.seed(1) 
# no s? si esto es necesario, pero creo que s? para que los datos de
# las simulaciones sean m?s reproducibles

################################
###### AVISO A NAVEGANTES ######
################################

## Para el primer caso no hay mucho problema. Las simulaciones corren sin tardar
## demasiado (de momento el ?nico par?metro que he elegido "razonando" es el 
## tiempo de simulaci?n, que es el que me dejaba el ordenador sin ir como una 
## patata) y los plots van bien, aunque algunos de simulaciones tardan un poco.

## Para el segundo ejemplo fuerzo un poco el l?mite de genotipos permitido en el
## evalAllGenotypes, por lo que simplifico un poco el DAG (est? explicado abajo)
## Las simulaciones corren y los plots salen, pero cuidado porque tiran un poco.

## Faltar?a ver si intento implementar las funciones
## oncoSimulPop/oncoSimulSample para hacer alguna cosa chula.

## Tama?o inicial - 2000, sigo papers Ram?n 2018 y 2015
## Tasa mutaci?n (mu) - 1e-5, seg?n Ram?n podr?a ser tambi?n 1e-6

## En cuanto a los plots, he hecho por n?mero de drivers, porque si no era una
## locura y sal?an ilegibles.

#############################################################
#### Ejemplo c?ncer colorrectal - datos Wood et al. 2007 ####
#############################################################

#### PathTiMEx encuentra dos posibles trayectorias para este conjunto
#### de mutaciones, dependiendo de si fuerza o no un "orden lineal"

#### Como son pocos genotipos, s? permite seguir adelante y hacer simulaciones


##  Caso I - Sin forzar progresi?n lineal del tumor
## Esquema de la Fig. 3A, panel izquierdo

## Lo interesante aqu? es que con OncoSimulR no podemos decir "no sabemos qu?
## requisitos tiene esta mutaci?n, pero sabemos que temporalmente ocurre despu?s
## de X". Esto es lo que dice realmente la Fig.3A con respecto a TP53/EVC2 y
## PIK3CA/EPHA3: no se conocen sus requisitos pero s? se especifica que ocurren 
## despu?s de KRAS y en estadios tard?os de la progresi?n tumoral respectivamente

## Opci?n 1 - nodos con requisitos desconocidos vienen directamente de la ra?z ##

## Como plotearemos por numero de drivers, creamos vector con sus nombres:
wood_drv <- c("APC", "KRAS", "TP53", "EVC2", "FBXW7", "TCF7L2", "PIK3CA", "EPHA3")

## Creacion del DAG de restricciones:
wood_root <- allFitnessEffects(data.frame(
                              parent = c("Root", "Root", "Root", "A", "B", "C"),
                              child = c("A", "C", "E", "B", "D", "D"),
                              s = 0.2,
                              sh = -0.9,
                              typeDep = "MN"),
                              geneToModule = c("Root" = "Root",
                                               "A" = "APC",
                                               "B" = "KRAS",
                                               "C" = "TP53, EVC2",
                                               "D" = "FBXW7, TCF7L2",
                                               "E" = "PIK3CA, EPHA3"),
                              drvNames = wood_drv
                                )
plot(wood_root, expandModules = TRUE)
eAG_wood_root <- evalAllGenotypes(wood_root)

plot(eAG_wood_root) #, use_ggrepel = TRUE) --> no funciona, demasiadas etiquetas
plot(eAG_wood_root, show_labels = FALSE) # no quita las etiquetas realmente??

# To try to see the genotypes that appear:
simul_wood_root <- oncoSimulIndiv(wood_root,
                              model="McFL",
                              onlyCancer = FALSE,
                              finalTime = 5000,
                              verbosity = 0,
                              mu = 1e-5,
                              initSize = 2000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              detectionProb = NA,
                              detectionSize = NA,
                              detectionDrivers = NA,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE
                              )

## Since plotting genotypes results in a very busy plot, we will only
## plot the number of drivers:

plot(simul_wood_root, addtot = TRUE, lwdClone = 0.9, log = "",
     thinData = TRUE, thinData.keep = 0.3,
     plotDiversity = TRUE,
     xlim = c(0, 5000)) 

plot(simul_wood_root, type = "stacked", thinData = TRUE, 
     thinData.keep = 0.1,
     plotDiversity = TRUE,
     xlim = c(0, 5000))

## We can see that the main clones have 5 mutations in driver genes, 
## which makes sense according to the DAG. Clones with 6-7 mutated drivers are
## those in which genes with mutual exclusivity are mutated. Their frequency is 
## very low, but they do appear: sh is not -1 because the most likely explanation
## for this mutual exclusivity is a lack of positive selection rather than 
## synthetic lethality (i.e. clones with both mutations are not viable.)

## He pedido que en los plots se muestre toda (de tiempo 0 a 5000) porque me 
## parece chulo ver lo que pasa con los clones con 5 vs. 6 drivers mutados

## If we look at the phylogenetic tree, we can see that there are clones in which
## the order APC --> KRAS --> TP53, which is well-established in the literature,
## is not respected (as we could already guess from the DAG of restrictions):

plotClonePhylog(simul_wood_root)


## Option 1B - Order effects ##

## We could use the same DAG of restrictions, but adding order effects
## The authors do not mention / suggest the existence of these effects in this
## case, but I believe this hypothesis is not contradictory to the graph either.
## In fact, I think the resulting graph is the one that makes the least assumptions
## about the data, even though it might not be the one that performs better in 
## cancer cell simulations (see discussion below).

## - The positive effect of mutating TP53/EVC2 before KRAS is very low. 
## - The same thing happens when mutating PIK3CA/EPHA3 before KRAS..
## - I do not penalize with a negative fitness effect because this would imply 
## that the clone is less viable when mutations happen this way, and I am not 
## sure that we would want to say that.

wood_order <- allFitnessEffects(data.frame(
        parent = c("Root", "A", "B", "C"),
        child = c("A", "B", "D", "D"),
        s = 0.2,
        sh = -0.9,
        typeDep = "MN"),
        geneToModule = c("Root" = "Root",
                         "A" = "APC",
                         "B" = "KRAS",
                         "C" = "TP53, EVC2",
                         "D" = "FBXW7, TCF7L2",
                         "E" = "PIK3CA, EPHA3"),
        drvNames = wood_drv,
        orderEffects = c("C > B" = 0.05,
                         "E > B" = 0.01) #, "E > C" = 0.01, "E > A" = 0.001),
)

plot(wood_order, expandModules = TRUE)

eAG_wood_order <- evalAllGenotypes(wood_order)

## c?digo simulaciones, por tenerlo por si acaso (no creo que sea especialmente
## interesante, yo me centrar?a en la discusi?n) (ver abajo)

# simul_wood_order <- oncoSimulIndiv(wood_order,
#                                    model="McFL",
#                                    onlyCancer = FALSE,
#                                    finalTime = 5000,
#                                    verbosity = 0,
#                                    mu = 1e-5,
#                                    initSize = 2000,
#                                    keepPhylog = TRUE,
#                                    seed = NULL,
#                                    detectionProb = NA,
#                                    detectionSize = NA,
#                                    detectionDrivers = NA,
#                                    errorHitMaxTries = FALSE,
#                                    errorHitWallTime = FALSE
# )
# 
# plotClonePhylog(simul_wood_order)
# 
# plot(simul_wood_order, addtot = TRUE, lwdClone = 0.9, log = "",
#      thinData = TRUE, thinData.keep = 0.3,
#      plotDiversity = TRUE,
#      xlim = c(0, 5000)) 
# 
# plot(simul_wood_order, type = "stacked", thinData = TRUE, 
#      thinData.keep = 0.1,
#      plotDiversity = TRUE,
#      xlim = c(0, 5000))

## Opci?n 2 - nodos con requisitos desconocidos vienen de KRAS

wood_kras <- allFitnessEffects(data.frame(
                              parent = c("Root", "A", "B", "B", "B", "C"),
                              child = c("A", "B", "D", "C", "E", "E"),
                              s = 0.2,
                              sh = -0.9,
                              typeDep = "MN"),
                              geneToModule = c("Root" = "Root",
                                               "A" = "APC",
                                               "B" = "KRAS",
                                               "C" = "TP53, EVC2",
                                               "D" = "PIK3CA, EPHA3",
                                               "E" = "FBXW7, TCF7L2"),
                              drvNames = wood_drv
                              )

plot(wood_kras, expandModules = TRUE)
eAG_wood_kras <- evalAllGenotypes(wood_kras)
plot(eAG_wood_kras)

# To try and see the genotypes that appear:
simul_wood_kras <- oncoSimulIndiv(wood_kras,
                                  model="McFL",
                                  onlyCancer = FALSE,
                                  finalTime = 5000,
                                  verbosity = 0,
                                  mu = 1e-5,
                                  initSize = 2000,
                                  keepPhylog = TRUE,
                                  seed = NULL,
                                  detectionProb = NA,
                                  detectionSize = NA,
                                  detectionDrivers = NA,
                                  errorHitMaxTries = FALSE,
                                  errorHitWallTime = FALSE)

## Again, the two possible plots:

plot(simul_wood_kras, addtot = TRUE, lwdClone = 0.9, log = "",
     thinData = TRUE, thinData.keep = 0.3,
     plotDiversity = TRUE,
     xlim = c(0, 5000)) 

plot(simul_wood_kras, type = "stacked", thinData = TRUE, 
     thinData.keep = 0.1,
     plotDiversity = TRUE,
     xlim = c(0, 1000))

# aqui no ganamos nada por poner hasta t = 5000 asi que lo he mantenido "corto"


##  Caso II - Forzando progresi?n lineal del tumor
## Esquema de la Fig. 4A, panel izquierdo

## This particular case was made by the authors to compare their approach to the
## one by Raphael and Vandin (2014). This model does not respect the order
## in driver mutations that we already mentioned and imposes an "artificial" 
## restriction on tumor progression (linear progression). However, we decided 
## to implement this example as well to see how it does in simulations.

## ^ Algo as? para decir que sabemos que en realidad esto no vale pa n?, pero 
## que ya puestos a modelar...

wood_linear <- allFitnessEffects(data.frame(parent = c("Root", "A", "B", "C"),
                                     child = c("A", "B", "C", "D"),
                                     s = 0.2,
                                     sh = -0.9,
                                     typeDep = "MN"),
                                  geneToModule = c("Root" = "Root",
                                                   "A" = "APC, EPHA3",
                                                   "B" = "EVC2, PIK3CA, TP53",
                                                   "C" = "KRAS, TCF7L2",
                                                   "D" = "FBXW7"),
                                  drvNames = wood_drv
                                  )


plot(wood_linear, expandModules = TRUE, autofit = TRUE)
eAG_wood_linear <- evalAllGenotypes(wood_linear)
plot(eAG_wood_linear) #, use_ggrepel = TRUE)

simul_wood_linear <- oncoSimulIndiv(wood_linear,
                              model="McFL",
                              onlyCancer = FALSE,
                              finalTime = 5000,
                              verbosity = 0,
                              mu = 1e-5,
                              initSize = 2000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              detectionProb = NA,
                              detectionSize = NA,
                              detectionDrivers = NA,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)


plot(simul_wood_linear, addtot = TRUE, lwdClone = 0.9, log = "",
     thinData = TRUE, thinData.keep = 0.3,
     plotDiversity = TRUE, 
     xlim = c(0, 1000))

plot(simul_wood_linear, type = "stacked", thinData = TRUE, 
     thinData.keep = 0.1,
     plotDiversity = TRUE,
     xlim = c(0, 1000))

## No s? si esto es relevante, pero lo dejo aunque sea para nosotros:

# cosa curiosa: con xlim = (0, 5000) salen poblaciones con 5 mutaciones --> 
# ?muestreo demasiado "largo"? (es decir, t demasiado grande?)

# Mirando el material suplementario del paper de Ram?n parece que s? estamos 
# haciendo un muestreo un poco largo si consideramos que lo que vemos es el 
# momento del diagn?stico. He decidido plotear la simulaci?n "completa" (hasta
# t = 5000 en xlim) solo en el primer ejemplo, para ense?ar lo que pasa, y
# en el resto mantener hasta t = 1000 (por ahorrar un poco de tiempo al plotear)

# En las figuras en las que hacemos hasta t = 5000 vemos que, a medida que crece
# el tumor, aparecen clones con m?s mutaciones, pero estas no parecen "fijarse"
# Yo creo que tiene sentido: el tumor sigue acumulando mutaciones, pero como 
# estas no aportan ninguna ventaja, estas poblaciones clonales no tienen ventaja
# y la presi?n selectiva no act?a a su favor (tampoco en su contra)

# --------------------- Discussion of the simulated data --------------------- # 

## We have seen four approaches. Three of them use Figure 3A as base to build
## the DAG, with small differences as there are two nodes (TP53/EVC2 and PIK3CA/
## EPHA2) that the authors describe as "independent" and that do not have a
## parent node in the graph. In order to include these nodes in our model, we 
## have devised three possible approaches: 

## (1) these nodes come from the root node (i.e. the wild type); 
## (2) these nodes are independent (as the authors suggest) and the effects of 
## these mutations are different depending on the order in which they occur;
## (3) these nodes are both child nodes from KRAS. 

## The last approach (4) is based on Figure 4A, where a linear progression is 
## imposed, but we will not focus on this one as we believe it is imposing 
## artificial restraints on tumor progression.

## We believe that approaches (2) and (3) are the ones that best model the 
## information provided by the authors, as they both consider the fact that 
## mutations in TP53/EVC2 and in PIK3CA/EPHA3 are observed at certain points in
## time. But which one is better? We selected some genotypes from the 
## evalAllGenotypes output in order to see exactly what is going on:

## Genotypes starting with PIK3CA (we could do this with more genotypes if
## we wanted to):
eAG_wood_root[grep("^PIK3CA", eAG_wood_root[,1]),]
# Genotype Fitness
# 6                PIK3CA   1.200
# 34       PIK3CA, TCF7L2   0.120
# 35         PIK3CA, TP53   1.440
# 92 PIK3CA, TCF7L2, TP53   0.144

eAG_wood_order[grep("^PIK3CA", eAG_wood_order[,1]),]
# Genotype Fitness
# 6                PIK3CA     1.0
# 34       PIK3CA, TCF7L2     0.1
# 35         PIK3CA, TP53     1.0
# 92 PIK3CA, TCF7L2, TP53     0.1

eAG_wood_kras[grep("^PIK3CA", eAG_wood_kras[,1]),]
# Genotype Fitness
# 6                PIK3CA   0.100
# 34       PIK3CA, TCF7L2   0.010
# 35         PIK3CA, TP53   0.010
# 92 PIK3CA, TCF7L2, TP53   0.001

## Which of these is better? Well, it seems that in wood_kras clones with the 
## genotype "PIK3CA" will appear less frequently, as their fitness is very low.
## Should the fitness of these clones be lower than 1? This is an interesting
## question. In this case, the fitness of 0.100 comes from the sh penalty applied
## to genotypes that deviate from the DAG of restrictions. 

## We probably need to know a lot more about the mechanisms of tumor progression
## in colorectal cancer to answer this question correctly. If we focus only on 
## replicating what the authors saw in this data set, the wood_kras model is 
## probably the "best" (most accurate) one, and the one that best reproduces
## the APC --> KRAS --> TP53 mutation order.

## Even when the authors could not find a dependency of having KRAS
## mutated in order to have either of the modules TP53/EVC2 or PIK3CA/EPHA2 
## mutated, which is what our DAG reflects, this simplification was necessary to 
## integrate the temporal information that the authors had, but could not 
## translate into a graph of dependencies. 


## Lo que quiero decir ah? arriba ^ es que s? que no es correcto decir que estos
## dos m?dulos dependen de KRAS porque no es eso lo que dice el paper, pero que 
## es la ?nica forma de integrar de alguna manera la informaci?n temporal que
## comentan (incluso as? perdemos algo de informaci?n).

## Creo que en realidad solo usan el modelo lineal para comparar su m?todo con 
## el de Raphael and Vandin, as? que si vemos que queda muy largo lo podemos 
## quitar de la presentaci?n (me parecer?a tonter?a quitarlo del documento que 
## entreguemos porque al final es trabajo hecho...)

###########################################
#### Ejemplo glioblastoma - Datos TCGA ####
###########################################

## Fig. 3C, panel izquierdo
## Volvemos a tener el problema del DAG vs. progresi?n temporal.

## Drivers: (las lineas estan x modulos mas o menos para editar facilmente)
gb_drv <- c("CDKN2A", "CDK4", "TP53-M", "MDM4", "MDM2", 
            "NF1", "FAF1", "SPTA1", "OBSCN", "CNTNAP2", 
            "PTEN-M", "PTEN-D", "PIK3CA", "IDH1",
            "PDGFRA", "LRP2", "EGFR", "RB1", "TP53-D", "PAOX")

## 1 - Modelando los m?dulos "hu?rfanos" como que vienen del Root:
TCGA_gb_root <- allFitnessEffects(data.frame(
                  parent = c("Root", "A", "Root", "Root", "B", "C", "D", "Root"),
                  child = c("A", "B", "C", "D", "E", "E", "E", "F"),
                  s = 0.2,
                  sh = -0.9,
                  typeDep = "OR"),
                  geneToModule = c("Root" = "Root",
                                   "A" = "CDKN2A, CDK4",
                                   "B" = "TP53-M, MDM4, MDM2",
                                   "C" = "NF1, FAF1, SPTA1, OBSCN, CNTNAP2",
                                   "D" = "PTEN-M, PTEN-D, PIK3CA, IDH1",
                                   "E" = "PDGFRA, LRP2",
                                   "F" = "EGFR, RB1, TP53-D, PAOX"),
                  drvNames = gb_drv)
plot(TCGA_gb_root, expandModules = TRUE, autofit = TRUE)


## No funciona evalAllGenotypes porque hay demasiados genotipos. Simplificamos:
## 1 - elimino modulo C porque los genes no vienen en el paper de Ciriello, y
## como vamos a enlazar con eso me parec?a que podr?a ser coherente
## 2 - en el resto de m?dulos, nos quedamos con los genes que vienen en Ciriello
## 3 - junto las dos mutaciones de PTEN (D) en una sola

## Con estos cambios parece que la cosa funciona sin demasiado problema
## borro el entorno de R para hacer hueco (no s? si hace falta):

rm(list = ls())
gb_drv_simple <- c("CDKN2A", "CDK4", "TP53-M", "MDM4", "MDM2", 
                    "PTEN", "PIK3CA",
                    "PDGFRA", "EGFR", "RB1")
TCGA_gb_root_simple <- allFitnessEffects(data.frame(
                  parent = c("Root", "A", "Root", "Root", "B", "D", "Root"),
                  child = c("A", "B", "D", "E", "E", "E", "F"),
                  s = 0.2,
                  sh = -0.9,
                  typeDep = "OR"),
                  geneToModule = c("Root" = "Root",
                                   "A" = "CDKN2A, CDK4",
                                   "B" = "TP53-M, MDM4, MDM2",
                                   "D" = "PTEN, PIK3CA",
                                   "E" = "PDGFRA",
                                   "F" = "EGFR, RB1, TP53-D"),
                  drvNames = gb_drv_simple)

plot(TCGA_gb_root_simple, expandModules = TRUE, autofit = TRUE)
eAG_TCGA_rs <- evalAllGenotypes(TCGA_gb_root_simple, max = 2100) 
# ^ lo hace y no tarda demasiado
plot(eAG_TCGA_rs) #horrible, podemos omitirlo yo creo

# mismos parametros que antes
simul_TCGA_rs <- oncoSimulIndiv(TCGA_gb_root_simple,
                                    model="McFL",
                                    onlyCancer = FALSE,
                                    finalTime = 5000,
                                    verbosity = 0,
                                    mu = 1e-5,
                                    initSize = 2000,
                                    keepPhylog = TRUE,
                                    seed = NULL,
                                    detectionProb = NA,
                                    detectionSize = NA,
                                    detectionDrivers = NA,
                                    errorHitMaxTries = FALSE,
                                    errorHitWallTime = FALSE)

# con tiempo = 1000 sale muy bonico, con 5000 ya cuesta un poco hacer gr?ficos
# y empiezan a aparecer "genotipos prohibidos" (6+ drivers mutados)

plot(simul_TCGA_rs, addtot = TRUE, lwdClone = 0.9, log = "",
     thinData = TRUE, thinData.keep = 0.3,
     plotDiversity = TRUE,
     xlim = c(0, 1000))
plot(simul_TCGA_rs, type = "stacked", thinData = TRUE, 
     thinData.keep = 0.1,
     plotDiversity = TRUE,
     xlim = c(0, 1000))

## 2 - Otro posible ?rbol:
TCGA_gb_2 <- allFitnessEffects(data.frame(
                    parent = c("Root", rep("A", 3), rep(c("B", "C", "D"), 2)),
                    child = c("A", "B", "C", "D", rep("E", 3), rep("F", 3)),
                    s = 0.2,
                    sh = -0.9,
                    typeDep = "OR"),
                    geneToModule = c("Root" = "Root",
                                     "A" = "CDKN2A, CDK4",
                                     "B" = "TP53-M, MDM4, MDM2",
                                     "C" = "NF1, FAF1, SPTA1, OBSCN, CNTNAP2",
                                     "D" = "PTEN-M, PTEN-D, PIK3CA, IDH1",
                                     "E" = "PDGFRA, LRP2",
                                     "F" = "EGFR, RB1, TP53-D, PAOX"))
plot(TCGA_gb_2, expandModules = TRUE, autofit = TRUE)

## Simplificamos para los siguientes pasos, con mismo criterio de antes:

TCGA_gb_2_simple <- allFitnessEffects(data.frame(
                    parent = c("Root", rep("A", 2), rep(c("B", "D"), 2)),
                    child = c("A", "B", "D", rep("E", 2), rep("F", 2)),
                    s = 0.2,
                    sh = -0.9,
                    typeDep = "OR"),
                    geneToModule = c("Root" = "Root",
                                     "A" = "CDKN2A, CDK4",
                                     "B" = "TP53-M, MDM4, MDM2",
                                     "D" = "PTEN, PIK3CA",
                                     "E" = "PDGFRA",
                                     "F" = "EGFR, RB1, TP53-D"),
                    drvNames = gb_drv_simple)

plot(TCGA_gb_2_simple, expandModules = TRUE, autofit = TRUE)
eAG_TCGA_2s <- evalAllGenotypes(TCGA_gb_2_simple, max = 2100) 
# ^ lo hace y no tarda demasiado
plot(eAG_TCGA_2s) #horrible, podemos omitirlo yo creo

# mismos parametros que antes
simul_TCGA_2s <- oncoSimulIndiv(TCGA_gb_2_simple,
                                model="McFL",
                                onlyCancer = FALSE,
                                finalTime = 5000,
                                verbosity = 0,
                                mu = 1e-5,
                                initSize = 2000,
                                keepPhylog = TRUE,
                                seed = NULL,
                                detectionProb = NA,
                                detectionSize = NA,
                                detectionDrivers = NA,
                                errorHitMaxTries = FALSE,
                                errorHitWallTime = FALSE)


plot(simul_TCGA_2s, addtot = TRUE, lwdClone = 0.9, log = "",
     thinData = TRUE, thinData.keep = 0.3,
     plotDiversity = TRUE,
     xlim = c(0, 1000))

plot(simul_TCGA_2s, type = "stacked", thinData = TRUE, 
     thinData.keep = 0.1,
     plotDiversity = TRUE,
     xlim = c(0, 1000))
