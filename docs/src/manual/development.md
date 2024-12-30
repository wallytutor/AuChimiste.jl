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
- Code is Julia, nothing else, no wrappers allowed.

By the above we mean that development must integrate simultaneous functionality development, documentation, and testing. Good formatting and readability are expected. Pull requests not meeting these criteria will be systematically refused.

Regarding the style, we make two major break-ups with standard Julia guidelines:

- First, functions and constant names have words separated by underscores, following Python-style. We hope in the future this more reasonable/human readable naming scheme will be accepted in Julia community. Other naming schemes such as the use of *Pascal-case* in structures are followed.

- The second point is that library functions are strictly declared through the use of `function` keyword, even for trivial functions. While *pure-functional* style of function declaration is fine in applications, being able to fold code and hide any implementation details when reviewing structure is an appreciated feature by the core developers.

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
