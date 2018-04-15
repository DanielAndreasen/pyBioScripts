# pyDNA
Tool to work with DNA sequences and see similarities.

## Create a DNA object
```python
from DNA import DNA
d1 = DNA('Name', 'AGTCTGNNN-NTGAAA-A')
d2 = DNA('DNA2', 'AGTTTGNNN-NTGAAMMA')
```

## Functionalities
```python
# Difference
res = d1 - d2

# Length work as expected:
print(len(d1))

# An "equalness" method to show how many percent are identical
print(d1.equalness(d2))

# Plot two sequences and show the difference too (requires matplotlib installed)
d1.plot(d2)
```
