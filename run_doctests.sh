# introduction is the only .rst file with unique doctest code
python -m doctest -v docs/source/introduction.rst

## test the other python files
for f in rhealpixdggs/*.py
do
  python -m doctest -v "$f"
done