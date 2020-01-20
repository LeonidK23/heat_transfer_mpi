to run the project local on computer:
```bash
mpirun -np 5 ./stencil
```
What do we need to compute the heat transfer on 2d grid of processors:
- slice blocks of the size (m, n) instead of (N, n)
- additionally exchange with neighbours on the top/ bottom and diagonal neighbours
- insert of received blocks (especially from the neighbours on the diagonal)
 
