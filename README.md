# moorea_filmies

Code repostitory to run analyses and generate figures and manuscript for Nitta et al. "Intergenerational Niche Differentiation in Filmy Ferns (Hymenophyllaceae)".

All code is in [R](https://cran.r-project.org/). The [targets package](https://wlandau.github.io/targets/index.html) is used to manage the workflow. 

To run all analyses and generate the manuscript:

1. [Clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
2. [Download and unzip the data](#data)
3. [Run `targets::tar_make()` in the provided docker container](#reproducible-analysis-with-docker)

## Data

For each of the links below, click on "Download Dataset" or "Download all", then place the zipped data file in the `data` folder of this repo. You can manually unzip the data archives if you want to see the contents, but the code needs the original zipped file in `data/` to run.

- https://figshare.com/s/d6349abf01a3756a5aae
- https://doi.org/10.5061/dryad.df59g
- https://doi.org/10.5061/dryad.fqz612jps

## Reproducible analysis with Docker

The analysis code requires various packages to be installed, and may not work properly if package versions have changed. Therefore, a 
[Docker image is provided](https://hub.docker.com/r/joelnitta/moorea_filmies) to run the code reproducibly. You can 
[install docker from here](https://docs.docker.com/install/).

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine):

```
cd /path/to/repo
```

Unzip the data:

```
docker run --rm -v ${PWD}:/tmpdir -w /tmpdir joelnitta/moorea_filmies:0.0.1 Rscript R/unzip_data.R
```

Run `targets::tar_make()`:

```
docker run --rm -v ${PWD}:/tmpdir -w /tmpdir joelnitta/moorea_filmies:0.0.1 Rscript -e 'targets::tar_make()'
```

You will see the targets being built by `targets`. The final manuscript should be compiled at the end as `manuscript.docx` (MS for journal submission) and `moorea_filmies_preprint.pdf` (preprint PDF) in the `results/ms` folder. Other figure and table files will also be compiled. Supplemental information will be written to the `results/si` folder.

## Interacting with the code

If you want to interact with the code in the Docker container, you can launch the container in the background using `docker-compose`:

```
docker-compose up -d
```

Navigate to http://localhost:8787/ in your browser of choice (firefox or google chrome recommended). There, you should be able to access an instance of the [RStudio](https://rstudio.com/) IDE, which can be used to inspect and manipulate objects in R.

When you're done, take down the container:

```
docker-compose down
```

## Licenses

- All code in this repository is licensed under the [MIT license](LICENSE)
- The [Roboto font](https://github.com/google/roboto/) is licensed under the [Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0)
