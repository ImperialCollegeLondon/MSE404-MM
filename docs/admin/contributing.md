## Installation guide

This website is hosted by Github pages and built by `mkdocs`. Under the hood it
uses github action to build the static page to `gh-pages` branch whenever a
new commit is found on the the `master` branch. This is done by
`.github/workflows/ci.yml`. In anycase if a deploy failed, please check the logs
of github action.

To contribute to the website, you need to install `mkdocs` and
`mkdocs-material` by following the instructions below:

Create and activate your virtual environment by:

```bash
 conda create -n mkdocs python=3.11
 conda activate mkdocs
```

Then install the following packages:
```bash
pip install mkdocs-material
pip install mkdocs-macros-plugin
```

Now you are ready to clone the repository:
```bash
git clone git@github.com:ImperialCollegeLondon/MSE404-MM.git
```
Note that you need to [:link:add your ssh-key to github](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)
before performing the clone.

Finally, go to the repo dir and spawn the server using:
```bash
mkdocs serve
```
then you should be able to see the compiled website at `http://127.0.0.1:8000`. 
Any change you make to the markdown files will be automatically updated in the
website so there's no need to re-do `mkdocs server`.

Now it's time to edit the files. The syntax is shown in the next section.

When you are done editing the files. Save the changes and deploy the website by:
```bash
# stage the change
git add .
# commit the change with message.
git commit -m "message."
```

Before pushing the changes, make sure that you have the latest version of the
`master` branch by:

```bash
git pull
```

Finally, push the changes to the remote repository by:

```bash
# push the change to the remote repository.
git push -u origin master
```
Note that you need to [:link:add your ssh-key to github](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)
before performing the push.

------------------------------------------------------------------------------

## Writing guide

The following information is a shorter version of the official documentation of
[:link:mkdocs-material](https://squidfunk.github.io/mkdocs-material/). You can
also find a similar one on our [:link:Group wiki](https://friendly-broccoli-22e4d939.pages.github.io/contributing/writing_wikis/#writing-wikis-basics).

- Standard Markdown Syntax:

```
# Title
## Subtitle
### Sub-subtitle
#### Sub-sub-subtitle

_italics_

**Bold**
```

- Adding a Hyperlink:
```
[something](assets/test.in)
```

- In-line code:
```
this is a `line`.
```

- Code blocks:
```
 ```python
 import numpy as np
 ```
```

- In-line math:
```
this is inline $\phi$ math
```

- Math blocks:
```
$$
\psi = \sum_a \phi_a
$$
```

- Numbered list:
```
1. test
2. test
```

    1. test
    2. test


- Dot list:
```
- something
- something
```

    - something
    - something

- Table:
```
| Method      | Description                          |
| ----------- | ------------------------------------ |
| `GET`       | :material-check:     Fetch resource  |
| `PUT`       | :material-check-all: Update resource |
| `DELETE`    | :material-close:     Delete resource |
```

    | Method      | Description                          |
    | ----------- | ------------------------------------ |
    | `GET`       | :material-check:     Fetch resource  |
    | `PUT`       | :material-check-all: Update resource |
    | `DELETE`    | :material-close:     Delete resource |

- Admonition (the indentation level N of admonition is determined by N*four spaces):
```
!!! note "show up"
    something.
	
??? note "hide, click to show up"
    something
```

    !!! note "show up"
        something.
    	
    ??? note "hide, click to show up"
        something
    

- Code annotation:
```
 ```python
 some code #(1)!
 ```

1.  annotation
```
    ```python
    some code #(1)!
    ```
   
    1.  annotation

- Code highlight (highlight frist and second line):
```
 ```python hl_lines="1-2"
 some code
 some code highlighted
 some other code highlighted
 some code
 ```
```
    ```python hl_lines="1-2"
    some code
    some code highlighted
    some other code highlighted
    some code
    ```

- Emojis (https://squidfunk.github.io/mkdocs-material/reference/icons-emojis/?h=emoji)
```
:link:
```
    :link:

- Images:
  - Center aligned (different syntax due to technical issues):
  ```
  <figure markdown="span">
    ![Diamond primitive cell](assets/Iron_bands.png)
  </figure>
  ```
      <figure markdown="span">
        ![Diamond primitive cell](https://dummyimage.com/600x400/)
      </figure>

  - Right aligned:
    ```
    ![MO](https://dummyimage.com/600x400/){ align=right}
    ```

    Some text. ![MO](https://dummyimage.com/600x400/){ align=right}



</br>

  - Resize:
    ```
    ![Image title](https://dummyimage.com/600x400/){ width="300" }
    ```


    ![Image title](https://dummyimage.com/600x400/){ width="300" }

- PNG files are prefered when showing structures of a material. You can export
  PNG with transparent background using VESTA.
