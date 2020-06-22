---
layout: page
title: Installation
permalink: /installation/
nav_order: 3
---

# Installation

Here we provide two ways to run our application. First is a standalone desktop application and the second is run the code directly on Rstudio

## Desktop Application

### Prerequisite

In order the application to run well in your Laptop/PC, you need have **JAVA environment** installed in your system. If you haven't install it yet, you can find it in [JAVA Website](https://www.java.com/). Meanwhile, you need have **Python environment** and **tensorflow package** installed in your system. If you haven't install it yet, we would recommend you to install Python through anaconda. You can find it in [Anaconda](https://www.anaconda.com), then install tensorflow package through by `pip install tensorflow` in terminal.

### For Mac

1. Download the `.dmg` file from Github.
2. Open the `.dmg` file and drag S2CNetDrug application to the desired location.
4. When opening the application for the first time, you may need to hold `^control`, click the application, click open, and go through verification for the application to open.

### For Windows

1. Download the `zip ` file from Github.
2. Unzip it and put it to the desired location.
3. Find `S2CNetDrug.exe` and click it to open application.

## Rstudio

If the latest version or more control on what is going on in application are desired, you can also download source code on Github and run it in Rstudio:

1. Download the source code file from Github.

2. Install all third-party packages in `denpendencies.R`. We provide code for you to download these packages. 

3. After all required packages are downloaded, open the `ui.R` or `server.R` and click `Run App` on the top right.