# Development

```@meta
CurrentModule = AuChimiste
```

This part of the documentation is intended for developers. It might also be useful for standard users trying to understand bugs or propose features. `AuChimiste` aims at having 100% first-level entities documented so that design features can be understood in the future.

## General guidelines

- Code written, code documented, code tested.
- Code lines make 72 characters, never more than 79.
- Code is not cluttered and comments are minimal.
- Code abuses of multiple dispatch if needed.
- Code is Julia, nothing else.

## Internals

Chemical Elements

```@docs
AuChimiste.ELEMENTS
AuChimiste.USER_ELEMENTS
AuChimiste.handle_element
AuChimiste.find_element
```

 Database parsing

```@docs
AuChimiste.DATA_PATH
AuChimiste.USER_PATH
```
