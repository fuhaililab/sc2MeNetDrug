---
layout: page
title: Installation
permalink: /installation/
nav_order: 3
---

# Installation

Here we provide three ways to run our application. First is a standalone desktop application and the second is run the code directly on Rstudio

## Docker

We also provide docker installation for sc2MeNetDrug. we recommend you to use Docker installation if possible. Since we set up all the environment appearpriately use docker image, so you do not worry about any environment issues that may happen during installation. For more information about docker installation, please send download request. 

## Desktop Application

### Prerequisite

In order the sc2MeNetDrug to run well in your Laptop/PC, you need have **JAVA environment** installed in your system. If you haven't install it yet, you can find it in [JAVA Website](https://www.java.com/). Meanwhile, you need have **Python environment** and **tensorflow package** installed in your system. If you haven't install it yet, we would recommend you to install Python through anaconda. You can find it in [Anaconda](https://www.anaconda.com), then install tensorflow package through by `pip install tensorflow` in terminal.

### For Mac

1. Download the `.dmg` file from Github.

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

## Rstudio

If the latest version or more control on what is going on in application are desired, you can also download source code on Github and run it in Rstudio:

1. Download the source code file from Github.

2. Install all third-party packages in `denpendencies.R`. We provide code for you to download these packages. 

3. After all required packages are downloaded, open the `ui.R` or `server.R` and click `Run App` on the top right.