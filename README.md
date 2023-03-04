# solar-system

This is an $N$-body simulation of the Solar System.

## Description

This is a first-persion view of the Solar System from Earth displaying the Moon, the Sun, and the
other seven planets. The Earth is looking directly at the Solar System's center of mass. This
direction forms the $\mathbf{\hat{z}}$ axis. The horizontal or $\mathbf{\hat{x}}$ axis is aligned so
that the span of $\mathbf{\hat{x}}$ and $\mathbf{\hat{z}}$ forms the ecliptic plane. The vertical or
$\mathbf{\hat{y}}$ axis is aligned so that the Earth orbits the Sun counterclockwise from the
perspective of an observer positioned on the positive $\mathbf{\hat{y}}$ axis looking towards the
origin.

This simulation isn't real-time. It's not even accurate. It's an $N$-body gravitational simulation.
The Sun starts out at the origin at rest. All of the planets and the Moon start at their approximate
real positions on 2023/01/01.

## Desktop and Web Browser

This application can be run on the desktop or in a web browser.

## Building It

The desktop version of solar-system is built using Cargo.

```console
cargo build --release
```

__TODO discuss building the web browser version__

## Future work

1. Document the code.
1. Finish creating unit tests.
1. Fix time keeping in WebAssemby version.
1. Fix label size issue in WebAssemby version.
1. Allow camera position and orientation to be changed through the UI.
1. Improve integration calculation of velocity and position when acceleration is changing rapidly.
1. Allow addition of celestial bodies through the UI.
1. Extend to support starting times other than 2023/01/01.
1. Replace measures library with a more mature library from crates.io.
1. Replace color library with conrod_core::color.
1. Fix bug in kiss3d that prevents text from being rendered in the correction location. This will be
   resolved if [kiss3d PR #319](https://github.com/sebcrozet/kiss3d/pull/319) is merged.
