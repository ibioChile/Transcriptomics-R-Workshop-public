# Workshop Setup

This workshop is meant to be hands-on and interactive. We would prefer not to spend time during the sessions getting everyone setup, so some things are needed in advance of the workshop. If you run into any problems, please feel free to email [Pamela Camejo](pcamejo@bio.puc.cl).

We will be using the software outlined in the instructions below. Please see the section for your operating system for those directions.

## Unix Shell

Unix Shell is a program that enables us to send commands to the computer and receive output. It is also referred to as the terminal or command line.

Computers with Mac OS or Linux already include a default Unix Shell program. The steps below describe some methods for identifying and opening a Unix Shell program if you already have one installed. There are also options for identifying and downloading a Unix Shell program, a Linux/UNIX emulator.

### Linux

The default Unix Shell for Linux OS is usually Bash. On most versions of Linux, it is accessible by running the (Gnome) Terminal or (KDE) Konsole or xterm, which can be found via the applications menu or the search bar.

### Mac OS

The default Unix Shell for a Mac computer is Bash. Your default shell is available via the Terminal program within your Utilities folder. To open Terminal, try one or both of the following:

- In Finder, select the Go menu, then select Utilities. Locate Terminal in the Utilities folder and open it.
- Use the Mac ‘Spotlight’ computer search function. Search for: Terminal and press Return.

Reference: [How to Use Terminal on a Mac](https://www.macworld.co.uk/how-to/mac-software/how-use-terminal-on-mac-3608274/)

### Windows

Computers with Windows operating systems do not automatically have a Unix Shell program installed. There are different solutions available for running Bash commands on Windows. 

- There is now a Bash shell command-line tool available for Windows 10, which we recommend to use in this workshop. This program is more than just a "terminal" since it can install **Linux within Windows in a safe way**. This allow us to locally run the majority of programs. This is known as Windows Subsystem for Linux. **Follow the instructions [in this link](https://docs.microsoft.com/en-us/windows/wsl/install-win10) for installation**.

- Additionally, you can run Bash commands on a remote computer or server that already has a Unix Shell, from your Windows machine. This can usually be done through a Secure Shell (SSH) client. One such client available for free for Windows computers is [PuTTY](https://www.putty.org/). Ask your PI (or senior bioinfo member of your lab group) if your lab has access to [the NLHPC](https://www.nlhpc.cl/). This is high-performance computing infraestructure free to use with research purposes. With the corresponding credentialas (username/password), you can access through a SSH client like Putty to  NLHPC linux servers and run all the commands of this workshop. 

- Another way is to use Virtual Machines. This options relies on windows native softwares that "inside" the software "emulate" a different machine (virtual machine). On this "virtual machine" we can run any Operative System like Ubuntu (Linux). You can find more info [into this link](https://www.howtogeek.com/170870/5-ways-to-run-linux-software-on-windows/). In case of Ubuntu, run the Xubuntu version that is a lightweight but fully functional

If none of the options below address your circumstances, try an online search for: Unix shell [your computer model] [your operating system].

# R

R and RStudio are separate downloads and installations. R is the underlying statistical computing environment and RStudio is a graphical integrated development environment (IDE) that makes using R much easier and more interactive. You need to install R before you install RStudio. Once installed, because RStudio is an IDE, RStudio will run R in the background. You do not need to run it separately.

### Windows

- Download R from the [CRAN website](http://cran.r-project.org/bin/windows/base/release.htm).
- Run the ```.exe``` file that was just downloaded
- Go to the [RStudio](https://www.rstudio.com/products/rstudio/download/#download) download page
- Under Installers select **RStudio x.yy.zzz - Windows Vista/7/8/10** (where x, y, and z represent version numbers)
- Double click the file to install it
- Once it’s installed, open RStudio to make sure it works and you don’t get any error messages.

### Mac OS

- Download R from the [CRAN website](http://cran.r-project.org/bin/macosx/).
- Select the ```.pkg``` file for the latest R version.
- Double click on the downloaded file to install R.
- Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download)
- Under Installers select **RStudio x.yy.zzz - Mac OS X 10.6+ (64-bit)** (where x, y, and z represent version numbers)
- Double click the file to install RStudio.
- Once it’s installed, open RStudio to make sure it works and you don’t get any error messages.

### Linux

- Follow the instructions for your distribution from [CRAN](https://cloud.r-project.org/bin/linux), they provide information to get the most recent version of R for common distributions. In any case, make sure you have at least R 3.2.
- Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download).
- Under Installers select the version that matches your distribution, and install it with your preferred method (e.g., with Debian/Ubuntu ```sudo dpkg -i rstudio-x.yy.zzz-amd64.deb``` at the terminal).
- Once it’s installed, open RStudio to make sure it works and you don’t get any error messages.

## Cytoscape

We will use Cytoscape to visualize and format transcriptional networks. Download it [here](https://cytoscape.org/download.html).
