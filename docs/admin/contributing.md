## Installation guide

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

Finally, go to the repo dir and spawn the server using:
```bash
mkdocs serve
```

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

- Dot list:
```
- something
- something
```

- Table:
```
| Method      | Description                          |
| ----------- | ------------------------------------ |
| `GET`       | :material-check:     Fetch resource  |
| `PUT`       | :material-check-all: Update resource |
| `DELETE`    | :material-close:     Delete resource |
```

- Admonition (the indentation level N of admonition is determined by N*four spaces):
```
!!! note "show up"
    something.
	
??? note "hide, click to show up"
    something
```

- Code annotation:
```
 ```python
 some code #(1)!
 ```

1.  annotation
```

- Code highlight (highlight frist and second line):
```
 ```python hl_lines="1-2"
 some code
 some code highlighted
 some other code highlighted
 some code
 ```
```

- Emojis (https://squidfunk.github.io/mkdocs-material/reference/icons-emojis/?h=emoji)
```
:link:
```

- Images:
  - Center aligned (different syntax due to technical issues):
  ```
  <figure markdown="span">
    ![Diamond primitive cell](assets/Iron_bands.png){ width="500" }
  </figure>
  ```

  - Right aligned:
    ```
    ![MO](assets/Mo_diagram.svg){: style="width:250px" align=right}
    ```
  - Resize:
    ```
    ![Image title](https://dummyimage.com/600x400/){ width="300" }
    ```

- PNG files are prefered when showing structures of a material. You can export
  PNG with transparent background using VESTA.
