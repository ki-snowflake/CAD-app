# Coronary Artery Disease (CAD) Translational Analytics Dashboard

An interactive Shiny dashboard exploring the link between clinical phenotypes and transcriptomic profiles in Coronary Artery Disease (CAD) research. Developed as a final project for the **Biostatistics** course.

🔗 **[View Live Application](https://snow8flake.shinyapps.io/cad-analysis-dashboard/)**

## 📌 Core Features

The dashboard includes 5 analytical tabs:

* **Patient Cohort & BMI:** Visualizes BMI distributions and cohort characteristics by clinical status.
* **Transcriptomic Profiling:** Evaluates differential expression of top 50 markers.
* **Functional Genomics (GO ORA):** Maps biological processes and gene-concept network linkages using `clusterProfiler`.
* **Comorbidity Risk Analysis:** Performs statistical assessments of risk factors, including Hypertension and Hyperlipidemia.
* **Biomarker Focus: CXCL5:** Statistically validates **CXCL5** as an atheroprotective factor, analyzing its correlation with disease progression.

<img width="2411" height="1507" alt="image" src="https://github.com/user-attachments/assets/16542c7b-99ae-416e-806a-2065e9157160" />

## 🧬 Key Research Insight

The analysis highlights the negative correlation between CXCL5 and CAD severity. The results suggest that CXCL5 acts as an independent endothelial survival factor, distinct from obesity surrogates (BMI), and that its downregulation serves as a marker for the exhaustion of vascular compensatory systems in advanced CAD stages.

## 🚀 Local Deployment

To run this application locally, clone the repository and execute the following in RStudio:

```bash
# Install required dependencies
install.packages(c("shiny", "tidyverse", "ggpubr", "viridis", "enrichplot", "shinythemes", "DT"))

# Launch the application
shiny::runApp()
```
