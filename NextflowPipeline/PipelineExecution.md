# MicrobioteRNApipe - Quality Control & Taxonomic Profiling Pipeline

Pipeline bioinformatique Nextflow pour le pré-traitement (Trimmomatic) et le profilage taxonomique (Kraken2 avec extraction des minimizers et kraken-biom) de données méta-transcriptomiques / RNA-Seq intratumorales.

---

## 1. Prérequis & Téléchargement des bases de données

### 1.1 Clonage du dépôt
```bash
git clone [https://github.com/GeNeHetX/MicrobioteRNApipe.git](https://github.com/GeNeHetX/MicrobioteRNApipe.git) 
cd MicrobioteRNApipe
```

### 1.2 Base de données Kraken2 PlusPFP
Téléchargez et décompressez l'archive de la base Kraken2 PlusPFP dans votre répertoire de références :

```bash
mkdir -p path_to/kraken_BD
cd path_to/kraken_BD

# Téléchargement de la base PlusPFP
wget [https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20260626.tar.gz](https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20260626.tar.gz)

# Décompression et nettoyage
tar -xvzf k2_pluspfp_20260626.tar.gz
rm k2_pluspfp_20260626.tar.gz
```
### 1.3 Fichier des adaptateurs Trimmomatic

Vérifiez la présence du fichier FASTA d'adaptateurs Illumina (utilisation du séquenceur Smarter mais à adapter selon le séquenceur):
``` bash
path_to/MicrobioteRNApipe/Nextflow_Pipeline/TruSeq3-PE-2-GGGGG.fa
```

# 2. Configuration (nextflow.config)

Le pipeline est configuré pour s'exécuter avec le gestionnaire SLURM et les conteneurs Singularity.

Les principaux paramètres définis dans la configuration :

- **krakendir** : Répertoire contenant les fichiers FASTQ d'entrée (UnmappedReads_output)

- **output_dir** : Répertoire de destination des résultats

- **kraken_db** : Chemin d'accès vers la base de données Kraken2 PlusPFP

- **confidence_score** : Seuil de score de confiance Kraken2 (ici 0.2)

- **adapter_path** : Fichier FASTA contenant les séquences d'adaptateurs pour Trimmomatic

# 3. Exécution du Pipeline Nextflow 

### 3.1 Exécution interactive via SLURM

Chargez les modules requis sur le cluster IFB et lancez l'analyse :

```bash
module load nextflow

nextflow run main_microbiote_pipeline.nf \
    -c nextflow.config \
    -profile slurm \
    -resume
```
### 3.2  Exécution en tâche de fond (Job SBATCH)

Pour exécuter le pipeline en tâche de fond sur le cluster IFB sans bloquer votre terminal, un script SLURM `launch_nf.job` est fourni.

**Structure du script** `launch_nf.job`
```bash
#!/bin/bash
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5GB
#SBATCH -c 1
#SBATCH -t 24:00:00

module purge
module load nextflow/26.04.4 

PIPELINE="/shared/projects/microbiote_pdacrna/MicrobioteRNApipe/NextflowPipeline/main_microbiote_pipeline.nf"
CONFIG="/shared/projects/microbiote_pdacrna/MicrobioteRNApipe/NextflowPipeline/test.config"

nextflow -C ${CONFIG} run${PIPELINE} -resume
```
**Lancement du script**

```bash
sbatch path_to/launch_nf.job
```

# 4. Extraction des métriques  (compter_reads_cat.sh)

Après l'exécution du workflow Nextflow, le script Bash `compter_reads_cat.sh` extrait les comptages de reads, les minimizers totaux et les minimizers uniques à partir de tous les fichiers `.kreport` générés, tout en assignant l'origine tissulaire de chaque échantillon.

L'assignation de l'origine se configure manuellement dans le script selon deux méthodes possibles :

---

### Option A : Définir un tableau d'assignation (`mapping`)
À utiliser si vos dossiers contiennent des échantillons issus de plusieurs origines tissulaires distinctes. Vous listez les noms de base (*basenames*) des échantillons associés à chaque modalité :

```bash
# Dans compter_reads_cat.sh
declare -A mapping

for id in 3_S17 11_S81 16_S26 17_S34; do mapping[$id]="Colon"; done
for id in 14_S10 27_S19 36_S91 45_S68; do mapping[$id]="Jejunum"; done
for id in R_CHC_1 R_CHC_2 R_CHC_3;   do mapping[$id]="CHC"; done
```

Voici comment extraire les basename par dossier pour créer la liste pour le for id ...:
```bash
ls path_to/kraken_calibration/score_0.2/*.kreport \
  | xargs -n 1 basename \
  | sed -E 's/^Unmapped_//; s/_[0-9]+\.[0-9]+\.kreport$//' \
  | sort -u \
  | tr '\n' ' '
```
### Option A bis : Générer directement la ligne complète
```bash
DOSSIER="path_to/kraken_calibration/score_0.2"
NOM_ORIGINE="ORIGINE"

IDS=$(ls "$DOSSIER"/*.kreport | xargs -n 1 basename | sed -E 's/^Unmapped_//; s/_[0-9]+\.[0-9]+\.kreport$//' | sort -u | tr '\n' ' ')

echo "for id in $IDS; do mapping[\$id]=\"$NOM_ORIGINE\"; done"
```

### Option B : Assigner une origine unique

À utiliser si l'ensemble des fichiers .kreport traités au cours de l'exécution provient de la même cohorte / du même tissu (ex. uniquement des échantillons IPMN). Il suffit de fixer la variable directement dans la boucle de traitement du script :

```bash
# Dans la fonction process_files() de compter_reads_cat.sh
origine="IPMN"
```
### 4.1 Rendre le script exécutable (dans BacterialEnrichmentAnalysis)
```bash
chmod +x compter_reads_cat.sh
```
### 4.2 Lancement sur les dossiers de résultats
Exécutez le script en passant en arguments un ou plusieurs répertoires contenant les fichiers .kreport :

```bash
./compter_reads_cat.sh dossier1 dossier2 dossier3
```
Exemple concret :
```bash
./compter_reads_cat.sh \
    /shared/projects/microbiote_pdacrna/AGASH/resultat_new_minimizer_kraken_TRIM_IPMN_RESTE/kraken_calibration/score_0.2 \
    /shared/projects/microbiote_pdacrna/AGASH/resultat_new_minimizer_kraken_TRIM_AUTRE_COHORTE/kraken_calibration/score_0.2
```
### 4.3 Fichier consolidé de sortie
Le script génère la table finale au format CSV, prête pour les analyses statistiques et la filtration sous R :

Fichier généré : /shared/projects/pipe_microbiote_pancreas/agash/test_extraction_TRIM_IPMN_TOUT.csv

Format des colonnes : Score,Sample,TaxID,Rank,Reads,Minimizers,Uniq_Minimizers,Categorie,ORIGINE,Methode,Species

### 5. Automatisation de la blacklist 

Avant toute analyse différentielle ou calcul de diversité, une liste d'exclusion (*blacklist*) est calculée dynamiquement via une **optimisation par courbes ROC (indice de Youden)**. Cette étape repose sur l'évaluation  de cohortes issues de tissus distincts pour distinguer le bruit technique du vrai signal biologique.



### 5.1 Cohortes requises (Multi-organes)

Pour garantir une calibration robuste, l'analyse doit intégrer des données issues de **localisation variés** préparés selon le même protocole :

* **Témoins digestifs / positifs :** `Colon` (sert de référence pour préserver le microbiote spécifique).
* **Organes / Tissus d'évaluation :** `CHC` (foie), `stomach_cancer` (estomac), `slice_pancreas_tumor_TT` (tranches pancréatiques), `pancreas_cancerPanNET` (tumeurs neuroendocrines), `ovary` (ovaire)



### 5.2 Principe de fonctionnement du script

1. **Annotation environnementale :** Associe les espèces bactériennes à leurs catégories d'isolation à partir de la base de données de référence (https://bacdive.dsmz.de/isolation-sources). Il faut télécharger la table entière en cliquant sur Download Table.
2. **Métriques combinées :**
   * **Prévalence inter-cohortes (`Presence_Controles_sans_IG`) :** Proportion de cohortes non digestives présentant le taxon.
   * **Couverture en minimizers (`Colon_S02_Moy_M`) :** Signal moyen de minimizers uniques dans le tissu de référence.
3. **Double seuillage ROC :**
   * Seuil de prévalence ubiquitaire maximal (Prévalence > Seuil_ROC).
   * Seuil de minimizers minimal (Minimizers < Seuil_ROC).


### 5.3 Lancement du script

Assurez-vous que les fichiers sources (`path_to/extraction_6_cohorte.csv` (créé par vos soins avec les 6 cohortes cités précédemment) et `fonctions_matrice.R`) sont accessibles, puis exécutez le script sous R :

```bash
source(path_to/blacklist.R)
```

### 5.4 Fichiers de sortie générés par la blacklist

L'exécution du script `blacklist.R` produit deux fichiers de référence essentiels dans votre répertoire de travail :

* **`Blacklist_strict.csv` :** Table des taxons classés comme contaminants/bruit systématique, contenant :
  * `Species` : Nom de l'espèce bactérienne à exclure.
  * `Groupe_Test` : Assignation de calibration (`Colon`, `Contaminant` ou `Autre`).
  * `Source_Microbiome` : Catégorie écologique issue de BacDive/DSMZ (`Microbiote Digestif/Oral`, `Contaminant Env`, `Autre Microbiote Humain/Cutané`, etc.).
  * `Presence_Controles_sans_IG` : Prévalence observée à travers les cohortes de contrôle non digestives.
  * `Colon_S02_Moy_M` : Moyenne de minimizers uniques par read dans la cohorte témoin.