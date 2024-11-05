# solar-system

This is an $N$-body simulation of the Solar System.

## Description

This is a first-person view of the Solar System from Earth displaying the Moon, the Sun, and the other seven planets. The camera is at the Earth looking directly at the Sun. This direction forms the $\mathbf{\hat{z}}$ axis. The horizontal or $\mathbf{\hat{x}}$ axis is aligned so that the span of $\mathbf{\hat{x}}$ and $\mathbf{\hat{z}}$ forms the ecliptic plane. The vertical or $\mathbf{\hat{y}}$ axis is aligned so that the Earth orbits the Sun counterclockwise from the perspective of an observer aligned with positive $\mathbf{\hat{y}}$ axis looking towards the origin.

This simulation isn't real-time. It's not even accurate. It's an $N$-body gravitational simulation. The Sun starts out at the origin at rest. All of the planets and the Moon start at their approximate real positions on 2023/01/01.

## Desktop and Web Browser

This application can be run on the desktop or in a web browser.

## Building It

The desktop version of solar-system is built using Cargo.

```console
cargo build --release
```

The Web Assembly version can be built using [Trunk](https://trunkrs.dev/).

```console
trunk build --release
```console

## Future work

1. Document the code.
1. Finish creating unit tests.
1. Fix time keeping in WebAssembly version.
1. Fix label size issue in WebAssembly version.
1. Allow camera position and orientation to be changed through the UI.
1. Improve integration calculation of velocity and position when acceleration is changing rapidly.
1. Allow addition of celestial bodies through the UI.
1. Extend to support starting times other than 2023/01/01.
