
pip install pyani-plus

pyani-plus --version

unzip cosmicus_genomes.zip

conda create -n ani -c bioconda -c conda-forge seaborn pyani-plus mummer

pyani-plus anim cosmicus_ani/ --database cosmicus.db --create-db --name cosmicus

pyani-plus export-run --database cosmicus.db --outdir cosmicus_ani/

pyani-plus plot-run --database cosmicus.db --outdir cosmicus_ani/