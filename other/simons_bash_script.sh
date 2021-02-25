#!/bin/bash

### .gitignore ###

# Create brand new .gitignore
rm .gitignore
touch .gitignore

# Fill .gitignore
echo "
# LaTeX
*.aux
*.fdb_latexmk
*.fls
*.log
*.out
*.synctex.gz

# Add every file larger than X MB to .gitignore
find -L . -size +1M | sed 's|^\./||g' | cat >> .gitignore
