---
layout: page
title: Installation
permalink: /installation/
nav_order: 3
---

# Installation

Here we provide three ways to run sc2MeNetDrug. Typically, sc2MeNetDrug will consume large memory during running if the size of scRNA-seq data is large. Regrad the normal size of scRNA-seq data, we suggest you to use the computer with at least 16g RAM(32g or more is recommended)

## Docker

We recommend you to use Docker installation if possible. Since we set up all the environment appearpriately use docker image, so you do not need to worry about any environment issues that may happen during installation and running. We provide the installation instruction in below, but please help us by filling in the form here [Download](./downloadRequest.md).

1. If you don't have docker environment in your computer, please download docker desktop first in [Docker website](https://www.docker.com). Make sure you set up docker environment appropriately.

2. Open the terminal/CMD, type following code to download sc2MeNetDrug image:

   `docker pull wfrain/sc2menetdrug`

3. After downloading, type following code to run the sc2MeNetDrug:

   `docker run -p 5000:5000 -v $HOME:$HOME wfrain/sc2menetdrug`

   Notice that we need to mount a computer directory by `-v` when using docker locally. Normally, In Mac, you can just run above code to mount local directory. If you using Windows, pleast check do you have `HOME` variable in path. Or, you can mount local directory manually using:

   `docker run -p 5000:5000 -v /C/Users:/Volumes wfrain/sc2menetdrug`

   This will mount local directory `/C/Users` to docker container directory `Volumes`. You can find it when you set working directory in sc2MeNetDrug.

   Meanwhile, To make sure sc2MeNetDrug run appropriately, please set the memory limit and swap to the maximum value in Docker desktop(Mac):

   <img src="../pic/docker.png" alt="docker" style="zoom:50%;" />

   Or you can set it when run the docker image:

   `docker run -p 5000:5000 -m 16g --memory-swap 4g -v /C/Users:/Volumes wfrain/sc2menetdrug`

   Where `-m` is the memory limit and `--memory-swap` is the amount of memory is allowed to swap to disk.

## Desktop Application

### Prerequisite

```
Important
{: .label .label-red }
```

In order the sc2MeNetDrug to run well in your Laptop/PC, you need have **JAVA environment** installed in your system. If you haven't install it yet, you can find it in [JAVA Website](https://www.java.com/). Meanwhile, you need have **Python environment** and **tensorflow package** installed in your **default Python environment**. If you haven't install it yet, we would recommend you to install Python through anaconda. You can find it in [Anaconda](https://www.anaconda.com), then install tensorflow package through by `pip install tensorflow` in terminal. The default Python environment refers to the Python environment that list in environment path. This is very important especially when your computer have more than one python environment.

### For Mac

1. Download the `.dmg` file from link.

2. Open the `.dmg` file and drag sc2MeNetDrug application to the desired location.

3. Before open the application, make sure you have download `gcc`. If not, you can download `gcc` by following steps: 

   1. Open the terminal, copy and run following code to install `homebrew`:

   `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"`

   2. Type and run following code to install `gcc`:`brew install gcc`.

4. When opening the application for the first time, you may need to hold `^control`, click the application, click open, and go through verification for the application to open.

### For Windows

1. Download the `zip ` file from Github.
2. Unzip it and put it to the desired location.
3. Find `sc2MeNetDrug.exe` and click it to open application.



## Video Demonstration

<iframe width="420" height="315" src="https://www.youtube.com/watch?v=q1BFcl_cdAk" frameborder="0" allowfullscreen></iframe>

