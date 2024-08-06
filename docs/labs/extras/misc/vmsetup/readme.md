Setting up the necessary software in a virtual machine
==============================================================

!!! note note
    If you are an MSE404 student you do not need to follow the instructions
    given here. 

These instructions are only for people who wish to run their own
copy of the software on their own machine. Doing is only recommended for more
advanced users or people without access to the compute server. It will mean many
files will not be at the location given in the lab documents so you will need to
be able to modify the instructions correctly yourself.

A virtual machine allows you to run a complete operating system contained
within an application regardless of what operating system you have. One
package that allows you to do this is [VMWare](https://www.vmware.com/).
VMWare is available free to Imperial College students via
<http://www.imperial.ac.uk/admin-services/ict/store/software/software-for-students/>.

Once you have it installed this, I suggest getting an Xubuntu iso from
<https://xubuntu.org/getxubuntu>, and using this to create a new virtual
machine. You can go with all the default settings, and pick whatever username
and password you like.

Once this is complete, it should reboot and prompt you for the password you
set during install to log you in to your desktop. Then first ensure the system
is fully up-to-date with

```bash
sudo apt-get update
sudo apt-get upgrade
```
You may also prefer to change the keyboard settings, since it defaults to a US
layout. If you go to the settings menu and pick "Language Support" you can
drag e.g. "English (United Kingdom)" to the top of the list and click the
"Apply System-Wide" button. You may also want to go to the "Regional Formats"
tab and select the UK option and apply system wide again. Then go to the
Keyboard settings menu and the "Layout" tab. At the bottom, click the "Add"
button and scroll down to UK English (there are several sub-options here that
may be better if you have e.g. a Mac). You can delete the US layout once you
have added the UK one.

Then you can simply install quantum espresso from the repositories, with

```bash
sudo apt-get install quantum-espresso
```

This will install a slightly older version of quantum espresso (which still
has all the features needed for these labs but may output additional files,
such as the wavefunction in a slightly different place than mentioned in the
lab description), and will not however create some of the directories, and
help files that we use in the examples. So some of the commands given in the
lab documents will not work without modification by you by default. In
particular, you should skip the various `module` commands as they will not be
necessary for you.

### If you have access to the server

One way around this, if you have access to the server, would be be to copy
these directories to the same place on your virtual machine. You can do this
by connecting to the college VPN and using the following commands,
substituting your college username for USER in the second command (which will
also prompt you for your college password):

```bash
sudo mkdir -p /opt/share
sudo rsync -auz USER@mt-studenty.mt.ic.ac.uk:/opt/share/quantum-espresso /opt/share
sudo mkdir /opt/Courses
sudo rsync -auvz USER@mt-studenty.mt.ic.ac.uk:/opt/Courses/MSE404 /opt/Courses
```

Then you should have everything necessary in place. I would suggest
re-running the last command at the beginning of each lab session to ensure
you have the latest version of the lab documents.

### If you do not have access to the server

There are two directories present on the server that are referred to at
various points in the lab documentation:

1. `/opt/share/quantum-espresso` - We have placed the various documentation
   files (txt, html and pdf) that are provided with the quantum-espresso
   source package in a directory here. It can be quite useful to refer to the
   text files from the command line. However the same information can all be
   found in the [documentation
   section](http://www.quantum-espresso.org/resources/users-manual) of the
   quantum espresso website. In particular the page with links to the [input
   data
   description](http://www.quantum-espresso.org/resources/users-manual/input-data-description).
2. `/opt/Courses/MSE404` - this lets students copy a directory for the lab
   containing the instructions and example files at the beginning of each lab
   session. You could clone the git repository
   [https://gitlab.com/eamonnmurray/MaterialsModelling.git](https://gitlab.com/eamonnmurray/MaterialsModelling.git)
   whereever you like and use this instead. If you are a student taking the
   course, I suggest pulling the repo at the beginning of each session as we
   make changes as the course proceeds each year.

  For example, you could clone the repository to the same directory as on
  the server:

```bash
sudo mkdir -p /opt/Courses
sudo git clone https://gitlab.com/eamonnmurray/MaterialsModelling.git /opt/Courses/MSE404
```

  This will create the folder `/opt/Courses/MSE404` with all the material
  for all the labs. To update this ahead of a lab session you can do

```bash
cd /opt/Courses/MSE404
sudo git pull
```
