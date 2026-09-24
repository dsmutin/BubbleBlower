# bubbleblower wiki

Iterative debubbling of totally coloured assembly graphs

**Status:** in development. Version: see repository file `VERSION` (starts at 0.0.1).

## Starting interconnections

```
user → bubbleblower CLI → run_pipeline() → JSON {status, ok, input_path}
                ↑
         tests (mandatory / optional)
                ↑
         examples/toy/run.py
```

Until real logic exists, `run_pipeline` is a baseline stub.

## Pages

- [Contracts](Contracts)
- [Testing](Testing)
