# vanni_et_al-figures
Code for the figures included in the manuscript [*Light into the darkness: Unifying the known and unknown coding sequence space in microbiome analyses*](https://www.biorxiv.org/content/10.1101/2020.06.30.180448v1)

To recreate the figures you can do as follow if you have a conda installation:

```bash
conda create -n vanni_et_al-figures
conda activate vanni_et_al-figures
conda install r=4.0 pkgconfig
```

Clone the repo with:

```bash
git clone https://github.com/functional-dark-side/vanni_et_al-figures.git
cd vanni_et_al-figures
```

Then get the data we used from [here](https://doi.org/10.6084/m9.figshare.12738476.v1). The data is stored in [SQlite](https://www.sqlite.org/index.html) dbs, one DB for each figure.

```bash
curl -JLO https://ndownloader.figshare.com/files/24109856
tar xvfz vanni_et_al-data.tar.gz
```

Then let's install the packages we used to plot the figures. First start R to get renv installed:

```
R
```

If everything went well, [renv](https://rstudio.github.io/renv/articles/renv.html) will be installed and you will get a message like:

```
* Installing renv 0.11.0 ... Done!
Successfully installed and loaded renv 0.11.0.
* Project '/vol/cloud/SANDBOX/vanni_et_al-figures' loaded. [renv 0.11.0]
```

And restore the environment:

```r
renv::restore()
q()
```

Now you can recreate the basic figures with:

```bash
./scripts/Figure1.R
```

> Figures will be saved in figures/

Or you can play with the code and recreate the figures in your favorite IDE. 

In the folder `figures-manuscipt` you can find the final figures used in the manuscript after some beautifying.
