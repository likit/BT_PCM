#!/usr/bin/env Rscript
Biological_Space <- read.csv("data/Biological_Space.csv", header=TRUE)
df <- data.frame(Biological_Space[,2:400])
df = df[, -nearZeroVar(df)]
set.seed(44321)
pca <- prcomp(df, retx=TRUE,scale.=TRUE)
scores <- pca$x[,1:3]
km <- kmeans(scores, center=3, nstart=5)
Protein <- Biological_Space[,1]
ggdata <- data.frame(scores, Cluster=km$cluster)
ggdata <- cbind(Protein, ggdata)
set.seed(32233)
pdf("Biological Space.pdf", width = 18, height = 18)
ggplot(ggdata) +
  geom_point(aes(x=PC1, y=PC2,
                 color=factor(Cluster)), size=5, shape=20) +
  ggtitle("Biological Space") +
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  geom_text(aes(x=PC1, y=PC2, label=Protein), size=2, hjust=0, alpha=0.6)
dev.off()
## Compounds
set.seed(42342)
Chemical_Space <- read.csv("data/Chemical_Space.csv", header=TRUE)
df <- data.frame(Chemical_Space[,2:45])
pca <- prcomp(df, retx=TRUE,scale.=TRUE)
set.seed(34233)
scores <- pca$x[,1:3]
km <- kmeans(scores, center=2, nstart=5)
Compound <- c("6,7-dimethoxy-1,3-dihydro-2-benzofuran-1-one",
              "(2S)-2-amino-3-phenylpropanoic acid",
              "(2S)-2-amino-3-(1H-imidazol-4-yl)propanoic acid",
              "1,2-diethyl benzene-1,2-dicarboxylate",
              "(prop-2-en-1-yl)thiourea",
              "5-nitro-1,3-thiazol-2-amine",
              "2,2,2-trichloroethane-1,1-diol",
              "(2R)-2-aminopropanoic acid",
              "3,4-dinitrobenzoic acid",
              "pentanedinitrile",
              "4-(pyridin-4-yl)pyridine",
              "pent-1-en-3-ol",
              "6-ethoxy-4,4-dimethyl-1,3-diazinane-2-thione",
              "2-(4-hydroxy-6-methyl-2-sulfanylidene-2,5-dihydropyrimidin-5-yl)acetate",
              "2-sulfanylidene-1,2,3,4-tetrahydropyrimidin-4-one",
              "2-{2-ethoxy-4-[(2,4,6-trioxo-1,3-diazinan-5-ylidene)methyl]phenoxy}acetic acid",
              "3-{3,5-dioxo-10-oxa-4-azatricyclo[5.2.1.0]dec-8-en-4-yl}benzoic acid",
              "[(E)-{amino[(4-fluorophenyl)amino]methylidene}amino]methanimidamide")
ggdata <- data.frame(scores, Cluster=km$cluster)
ggdata <- cbind(Compound, ggdata)
set.seed(10093)
pdf("Chemical Space.pdf", width = 14, height = 14)
ggplot(ggdata) +
  geom_point(aes(x=PC1, y=PC2,
                 color=factor(Cluster)), size=5, shape=20) +
  ggtitle("Chemical Space") +
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  geom_text(aes(x=PC1, y=PC2, label=Compound), size=2, hjust=0, alpha=0.6)
dev.off()
