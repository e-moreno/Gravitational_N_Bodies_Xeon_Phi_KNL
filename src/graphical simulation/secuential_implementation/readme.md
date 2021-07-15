# Compile

``` bash
    gcc -o nbodies_plot.exe nbodies_final.c -lm
```

# Generate data

``` bash
    ./nbodies_plot.exe 256 1024 > plot_data.csv
```

# Plot data

``` bash
    ./plot.py plot_data.csv
```