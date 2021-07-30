## 2D Ising Model in C with multi-threading

A multi-threading implementation of the Ising model on a 2D grid in the C programming language with as goal to make it as computational efficient as possible. In benchmarks i'm able to get 250+ million iteration on 6 threads. 

Video with increasing temperature: https://www.youtube.com/watch?v=_rn17hP78fY

Meme Video with Bad Apple!! in Ising: https://www.youtube.com/watch?v=tlGv9jt3gJg

## Usage

Change the definition in `ising-base.c` to the desired situation and run the following

```
gcc ising-base.c -lpthread
./a.exe
```

This will create a directory with frames of the 2D Ising grid. This is than analysed using the tools in `mypythontools.py` which are used in `base-analysis.py`. 