# bamread

Read bam files quickly into dataframes in Python

## Development

### Development

See [How to Publish an Open-Source Python Package to PyPI](https://realpython.com/pypi-publish-python-package/)

Install

```bash
pip install -e .  # -e flag means editing affects the installed package
```

Create sdist and wheel (requires [build](https://pypa-build.readthedocs.io/en/stable/))

```bash
python -m build
```

Check whether the package archives look fine:

```bash
twine check dist/*
```

Upload to PyPI test server

```bash
twine upload -r testpypi dist/*
```

Chek that the package can be installed from the test server

```bash
python -m pip install -i https://test.pypi.org/simple bamread
python -c "import bamread; print(bamread.read_bam('Some.bam'))"
```

Upload to pypi

```bash
twine upload -r testpypi dist/*
```

# Run github actions locally

See [act](https://github.com/nektos/act)

```bash
act -l # list available jobs
act -j <job> # run a job
act pull_request # run jobs triggered upon pull_request
```
