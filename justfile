@default:
    just --list --unsorted


# Aliases/shorthands for the recipes below
alias rs := rust
# alias readme := make-readme
alias zip := compress_html
alias doc := docs
alias qmd := quarto
alias env := setup-env
alias local := setup-env
alias dev := setup-env
alias all_docker := docker
alias container_prep := docker
alias py := python
alias all := preflight
alias doit := preflight


# Rust recipes
# ----------------------------------------------------

export CARGO_TERM_COLOR := "always"

# make sure rust is installed and show the toolchain version
rustup-show:
  rustup show

# run all Cargo checks
check:
  cargo check --all-targets

# Format all Rust files strictly
fmt:
  cargo fmt --all -- --check

# Run through all cargo clippy lints
clippy:
  cargo clippy --all-targets

# Run all available tests throughout the workspace
test:
  cargo insta test --workspace

# install a couple cargo extensions
install:
  cargo install --locked cargo-sort cargo-audit

# run Cargo sort on the whole workspace
sort: install
  cargo sort --check --workspace

# run cargo audit on the whole workspace
audit: install
  cargo audit

# Run all Rust recipes
rust: check fmt clippy test sort audit


# Quarto recipes
# ----------------------------------------------------

# render quarto documents
render:
    quarto render docs/index.qmd
    quarto render docs/examples.qmd

# turn rendered index quarto file into the repo readme
# make-readme:
# @mv docs/index.md ./README.md

# render developer guide
render-dev:
    quarto render docs/developer.qmd

# compress rendered HTML
compress_html:
    @gzip -f docs/index.html
    @gzip -f docs/developer.html
    @gzip -f docs/examples.html

quarto: render render-dev

docs: render render-dev compress_html # make-readme


# Docker recipes
# ----------------------------------------------------

# build docker image
docker-build:
    docker build .

# push docker image to Docker Hub
docker-push:
    docker push

docker: docker-build docker-push


# Python/full environment recipes
# ----------------------------------------------------

# install uv and use it to install PyPI dependencies
setup-env:
    curl -LsSf https://astral.sh/uv/install.sh | sh
    uv venv
    source activate .venv/bin/activate
    uv pip install -r requirements.txt

# install pre-commit hooks and run them
pre-commit: setup-env
    source activate .venv/bin/activate
    pre-commit install
    pre-commit run

# freeze current dependencies into an explicit lock file
py-freeze: setup-env
    source activate .venv/bin/activate
    uv pip freeze > requirements.txt

# run lints with Ruff
py-lints: setup-env
    source activate .venv/bin/activate
    ruff check . --exit-zero --fix --unsafe-fixes

# run Python formatting with Ruff
py-format: setup-env
    source activate .venv/bin/activate
    ruff format .

# Sort imports alphabetically
py-sort-imports: setup-env
    source activate .venv/bin/activate
    ruff check . -n --select=I --fix

python: py-lints py-format py-sort-imports py-freeze pre-commit


# Full pipeline recipes
# ----------------------------------------------------

format: fmt sort py-format py-sort-imports
lint: clippy py-lints
preflight: rust quarto docker python
