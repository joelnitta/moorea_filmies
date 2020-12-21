# moorea_filmies

Code repostitory to run analyses and generate figures and manuscript for Nitta et al. "Intergenerational Niche Differentiation in Filmy Ferns (Hymenophyllaceae)".

All code is in [R](https://cran.r-project.org/). The [targets package](https://wlandau.github.io/targets/index.html) is used to manage the workflow. 

To run all analyses and generate the manuscript:

1. [Clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
2. [Download and unzip the data](#data)
3. [Run `targets::tar_make()` in the provided docker container](#reproducible-analysis-with-docker)

## Data

For each of the links below, click on "Download Dataset", then place the zipped data file (it will have a similar name to the DOI) in the `data` folder of the 
this repo. You can manually unzip the data archives if you want to see the contents, but the code needs the original zipped file in `data/` to run.

- https://doi.org/10.5061/dryad.df59g
- https://doi.org/10.5061/dryad.fqz612jps

## Reproducible analysis with Docker

The analysis code requires various packages to be installed, and may not work properly if package versions have changed. Therefore, a 
[Docker image is provided](https://hub.docker.com/r/joelnitta/moorea_filmies) to run the code reproducibly. You can 
[install docker from here](https://docs.docker.com/install/).

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine), and launch the container:

```
cd /path/to/repo
docker-compose up -d
```

Run `targets::tar_make()` inside the container:

```
docker exec moorea_filmies_analysis_1 Rscript -e "targets::tar_make()"
```

You will see the targets being built by `targets`, and the final manuscript should be compiled at the end as `manuscript.pdf` and `manuscript.docx` in the `results/ms` folder. Other figure and table files will also be compiled.

When it's finished, take down the container:

```
docker-compose down
```

## Licenses

- All code in this repository is licensed under the [MIT license](LICENSE)
- The [Roboto font](https://github.com/google/roboto/) is licensed under the [Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0)
