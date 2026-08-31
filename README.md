# MicrobioteRNApipe - Pipeline de décontamination et de profilage taxonomique pour RNA-Seq FFPE à faible biomasse

## Description

**MicrobioteRNApipe** est un pipeline bioinformatique développé avec **Nextflow** et **Docker** pour le contrôle qualité, l'élimination des séquences de l'hôte, le profilage taxonomique et la décontamination calibrée par courbes ROC des données de microbiome intratumoral (notamment sur tissus FFPE à faible biomasse)

Le workflow traite les reads d'ARN non humains selon les étapes suivantes :
1. **Filtration de l'hôte :** Alignement contre le génome humain de référence via **STAR** (utilisation des lectures non alignées *unmapped* pour l'analyse microbienne)
2. **Prétraitement & Nettoyage :** Suppression des adaptateurs et filtration des bases de faible qualité via **Trimmomatic**
3. **Classification taxonomique :** Assignation ultra-rapide par $k$-mers et *minimizers* via **Kraken2** (avec la base de données exhaustive **PlusPFP**)
4. **Extraction des métriques & Consolidation :** Script utilitaire Bash (`compter_reads_cat.sh`) extrayant le nombre de reads, les minimizers totaux et les minimizers uniques par échantillon.
5. **Décontamination calibrée par ROC :** Élimination systématique des effets batch, du *kitome* et des erreurs d'assignation algorithmiques (`blacklist.R`)

---

## Auteur & Encadrement

- **Agash UTHAYAKUMAR** ([GitHub](https://github.com/AGASH0))
- Encadré par : **Camille Pignolet** et **Lucie Gomes** (Centre de Recherche sur l'Inflammation, Inserm)

---

## Dépendances & Environnement

Le pipeline est préconfiguré pour s'exécuter via le gestionnaire de jobs **SLURM** et les conteneurs **Singularity** sur cluster de calcul HPC (ex. cluster IFB Core) :

* [Nextflow](https://www.nextflow.io/) (23.10)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (v0.39)
* [STAR](https://github.com/alexdobin/STAR) (v2.7.3)
* [Kraken2](https://github.com/DerrickWood/kraken2) (v2.17.1 avec option `--report-minimizer-data`)
* [R](https://www.r-project.org/) (4.2) avec les packages `tidyverse`, `pROC`, `moments` et `MaAsLin2`

---

## Bases de données requises

Avant l'exécution, assurez-vous de disposer des bases suivantes sur votre espace de stockage :

1. **Index génomiques & transcriptomiques humains :** Préparés pour la filtration des reads de l'hôte avec STAR
2. **Base de données Kraken2 PlusPFP :** Base de référence incluant Bactéries, Archées, Virus, Protozoaires, Champignons, Plantes et Humain (Téléchargement : [k2_pluspfp_20260626.tar.gz](https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20260626.tar.gz))
3. **Table d'isolation des souches BacDive / DSMZ :** Métadonnées écologiques pour l'annotation des catégories d'habitat bactérien (Disponible sur [BacDive](https://bacdive.dsmz.de/isolation-sources)).

---

## Guide de démarrage rapide

### 1. Configuration des paramètres (`nextflow.config`)
Éditez votre fichier de configuration pour indiquer vos chemins de fichiers d'entrée, de sortie et de bases de données :

```groovy
params {
    krakendir         = "/chemin/vers/UnmappedReads_output"
    output_dir        = "/chemin/vers/resultats"
    kraken_db         = "/chemin/vers/kraken_BD"
    confidence_score  = 0.2
    adapter_path      = "/chemin/vers/TruSeq3-PE-2-GGGGG.fa"
}
```

### 2. Soumission de l'exécution Nextflow

Sur un cluster HPC géré par SLURM, lancez le pipeline en arrière-plan à l'aide du script fourni[cite: 3] :

```bash
sbatch launch_nf.job
```

### 3. Extraction des métriques & Décontamination ROC
Une fois le workflow terminé, consolidez les rapports .kreport et appliquez la calibration de la blacklist :

```bash
# Visualisation des comptages et des minimizers uniques
./compter_reads_cat.sh /chemin/vers/kraken_calibration/score_0.2

# Génération de la blacklist par optimisation ROC sous R
Rscript blacklist.R
```

# Documentation détaillée

Pour le guide complet étape par étape (configuration des tableaux de mapping, options d'exécution et tutoriel pas-à-pas), consultez le fichier PipelineExecution.md.


