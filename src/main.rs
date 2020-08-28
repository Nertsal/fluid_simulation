#![allow(dead_code)]

use piston_window::types::Color;
use piston_window::*;

const SCALE: f64 = 10.0;

fn main() {
    let n: usize = 64;

    let mut window: PistonWindow = WindowSettings::new(
        "Fluid simulation",
        [to_coord_u32(n as i32), to_coord_u32(n as i32)],
    )
    .exit_on_esc(true)
    .build()
    .unwrap();

    let mut fluid = Fluid::new(n, 16, 0.0000000001, 0.0001);

    while let Some(event) = window.next() {
        if let Some([mouse_x, mouse_y]) = event.mouse_cursor_args() {
            let fluid_x = (mouse_x / SCALE) as usize;
            let fluid_y = (mouse_y / SCALE) as usize;
            fluid.add_dye(fluid_x, fluid_y, 0.05);
            fluid.add_velocity(fluid_x, fluid_y, 0.0, 1.0)
        }
        window.draw_2d(&event, |c, g, _| {
            clear([0.0, 0.0, 0.0, 0.0] as Color, g);
            fluid.draw(&c, g);
        });

        event.update(|arg| {
            fluid.update(arg.dt);
        });
    }
}

fn to_coord(game_coord: i32) -> f64 {
    (game_coord as f64) * SCALE
}

fn to_coord_u32(game_coord: i32) -> u32 {
    to_coord(game_coord) as u32
}

pub fn draw_rectangle(
    color: Color,
    x: i32,
    y: i32,
    width: i32,
    height: i32,
    con: &Context,
    g: &mut G2d,
) {
    let gui_x = to_coord(x);
    let gui_y = to_coord(y);

    rectangle(
        color,
        [
            gui_x,
            gui_y,
            SCALE * (width as f64),
            SCALE * (height as f64),
        ],
        con.transform,
        g,
    );
}

struct Fluid {
    n: usize,
    iter: usize,

    diff: f64,
    visc: f64,

    s: Vec<f64>,
    density: Vec<f64>,

    vx: Vec<f64>,
    vy: Vec<f64>,

    vx0: Vec<f64>,
    vy0: Vec<f64>,
}

impl Fluid {
    fn new(n: usize, iter: usize, diff: f64, visc: f64) -> Self {
        Self {
            n,
            iter,
            diff,
            visc,
            s: vec![0.0; n * n],
            density: vec![0.0; n * n],
            vx: vec![0.0; n * n],
            vy: vec![0.0; n * n],
            vx0: vec![0.0; n * n],
            vy0: vec![0.0; n * n],
        }
    }

    fn add_dye(&mut self, x: usize, y: usize, amount: f64) {
        let index = Self::coords_to_index(self.n, x, y);
        self.density[index] += amount;
    }

    fn add_velocity(&mut self, x: usize, y: usize, amount_x: f64, amount_y: f64) {
        let index = Self::coords_to_index(self.n, x, y);
        self.vx[index] += amount_x;
        self.vy[index] += amount_y;
    }

    fn update(&mut self, delta_time: f64) {
        Self::diffuse(self.n, self.iter, 1, &mut self.vx0, &self.vx, self.visc, delta_time);
        Self::diffuse(self.n, self.iter, 2, &mut self.vy0, &self.vy, self.visc, delta_time);

        Self::project(
            self.n,
            self.iter,
            &mut self.vx0,
            &mut self.vy0,
            &mut self.vx,
            &mut self.vy,
        );

        Self::advect(
            self.n,
            1,
            &mut self.vx,
            &self.vx0,
            &self.vx0,
            &self.vy0,
            delta_time,
        );
        Self::advect(
            self.n,
            2,
            &mut self.vy,
            &self.vy0,
            &self.vx0,
            &self.vy0,
            delta_time,
        );

        Self::project(
            self.n,
            self.iter,
            &mut self.vx,
            &mut self.vy,
            &mut self.vx0,
            &mut self.vy0,
        );

        Self::diffuse(
            self.n,
            self.iter,
            0,
            &mut self.s,
            &self.density,
            self.diff,
            delta_time,
        );
        Self::advect(
            self.n,
            0,
            &mut self.density,
            &self.s,
            &self.vx,
            &self.vy,
            delta_time,
        );
    }

    fn draw(&self, con: &Context, g: &mut G2d) {
        for i in 0..self.n {
            for j in 0..self.n {
                let dye_amount = self.density[Self::coords_to_index(self.n, i, j)] as f32;
                let color: Color = [dye_amount, dye_amount, dye_amount, 1.0];
                draw_rectangle(color, i as i32, j as i32, 1, 1, con, g);
            }
        }
    }

    fn coords_to_index(n: usize, x: usize, y: usize) -> usize {
        x + y * n
    }

    fn diffuse(
        n: usize,
        iter: usize,
        b: u32,
        x: &mut [f64],
        x0: &[f64],
        diff: f64,
        delta_time: f64,
    ) {
        let a = delta_time * diff * ((n as f64) - 2.0) * ((n as f64) - 2.0);
        Self::lin_solve(n, iter, b, x, x0, a, 1.0 + 6.0 * a);
    }

    fn lin_solve(n: usize, iter: usize, b: u32, x: &mut [f64], x0: &[f64], a: f64, c: f64) {
        let c_recip = 1.0 / c;
        for _ in 0..iter {
            for j in 1..n - 1 {
                for i in 1..n - 1 {
                    x[Self::coords_to_index(n, i, j)] = (x0[Self::coords_to_index(n, i, j)]
                        + a * (x[Self::coords_to_index(n, i + 1, j)]
                            + x[Self::coords_to_index(n, i - 1, j)]
                            + x[Self::coords_to_index(n, i, j + 1)]
                            + x[Self::coords_to_index(n, i, j - 1)]))
                        * c_recip;
                }
            }
        }
        Self::set_bnd(n, b, x);
    }

    fn project(
        n: usize,
        iter: usize,
        vel_x: &mut [f64],
        vel_y: &mut [f64],
        p: &mut [f64],
        div: &mut [f64],
    ) {
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                div[Self::coords_to_index(n, i, j)] = (-0.5
                    * (vel_x[Self::coords_to_index(n, i + 1, j)]
                        - vel_x[Self::coords_to_index(n, i - 1, j)]
                        + vel_y[Self::coords_to_index(n, i, j + 1)]
                        - vel_y[Self::coords_to_index(n, i, j - 1)]))
                    / (n as f64);
                p[Self::coords_to_index(n, i, j)] = 0.0;
            }
        }
        Self::set_bnd(n, 0, div);
        Self::set_bnd(n, 0, p);
        Self::lin_solve(n, iter, 0, p, div, 1.0, 6.0);

        for j in 1..n - 1 {
            for i in 1..n - 1 {
                vel_x[Self::coords_to_index(n, i, j)] -= 0.5
                    * (p[Self::coords_to_index(n, i + 1, j)]
                        - p[Self::coords_to_index(n, i - 1, j)])
                    * (n as f64);
                vel_y[Self::coords_to_index(n, i, j)] -= 0.5
                    * (p[Self::coords_to_index(n, i, j + 1)]
                        - p[Self::coords_to_index(n, i, j - 1)])
                    * (n as f64);
            }
        }
        Self::set_bnd(n, 1, vel_x);
        Self::set_bnd(n, 2, vel_y);
    }

    fn advect(
        n: usize,
        b: u32,
        d: &mut [f64],
        d0: &[f64],
        vel_x: &[f64],
        vel_y: &[f64],
        delta_time: f64,
    ) {
        let dt_x = delta_time * ((n as f64) - 2.0);
        let dt_y = delta_time * ((n as f64) - 2.0);

        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let tmp1 = dt_x * vel_x[Self::coords_to_index(n, i, j)];
                let tmp2 = dt_y * vel_y[Self::coords_to_index(n, i, j)];
                let mut x = (i as f64) - tmp1;
                let mut y = (j as f64) - tmp2;

                if x < 0.5 {
                    x = 0.5
                }
                if x > (n as f64) + 0.5 {
                    x = (n as f64) + 0.5
                }
                let i0 = x.floor();
                let i1 = i0 + 1.0;
                if y < 0.5 {
                    y = 0.5
                }
                if y > (n as f64) + 0.5 {
                    y = (n as f64) + 0.5
                }
                let j0 = y.floor();
                let j1 = j0 + 1.0;

                let s1 = x - i0;
                let s0 = 1.0 - s1;
                let t1 = y - j0;
                let t0 = 1.0 - t1;

                d[Self::coords_to_index(n, i, j)] = s0
                    * (t0 * d0[Self::coords_to_index(n, i0 as usize, j0 as usize)]
                        + t1 * d0[Self::coords_to_index(n, i0 as usize, j1 as usize)])
                    + s1 * (t0 * d0[Self::coords_to_index(n, i1 as usize, j0 as usize)]
                        + t1 * d0[Self::coords_to_index(n, i1 as usize, j1 as usize)]);
            }
        }
        Self::set_bnd(n, b, d);
    }

    fn set_bnd(n: usize, b: u32, x: &mut [f64]) {
        for i in 1..n - 1 {
            x[Self::coords_to_index(n, i, 0)] = if b == 2 {
                -x[Self::coords_to_index(n, i, 1)]
            } else {
                x[Self::coords_to_index(n, i, 1)]
            };
            x[Self::coords_to_index(n, i, n - 1)] = if b == 2 {
                -x[Self::coords_to_index(n, i, n - 2)]
            } else {
                x[Self::coords_to_index(n, i, n - 2)]
            };
        }
        for j in 1..n - 1 {
            x[Self::coords_to_index(n, 0, j)] = if b == 1 {
                -x[Self::coords_to_index(n, 1, j)]
            } else {
                x[Self::coords_to_index(n, 1, j)]
            };
            x[Self::coords_to_index(n, n - 1, j)] = if b == 1 {
                -x[Self::coords_to_index(n, n - 2, j)]
            } else {
                x[Self::coords_to_index(n, n - 2, j)]
            };
        }

        x[Self::coords_to_index(n, 0, 0)] =
            0.5 * (x[Self::coords_to_index(n, 1, 0)] + x[Self::coords_to_index(n, 0, 1)]);
        x[Self::coords_to_index(n, 0, n - 1)] =
            0.5 * (x[Self::coords_to_index(n, 1, n - 1)] + x[Self::coords_to_index(n, 0, n - 2)]);
        x[Self::coords_to_index(n, n - 1, 0)] =
            0.5 * (x[Self::coords_to_index(n, n - 2, 0)] + x[Self::coords_to_index(n, n - 1, 1)]);
        x[Self::coords_to_index(n, n - 1, n - 1)] = 0.5
            * (x[Self::coords_to_index(n, n - 2, n - 1)]
                + x[Self::coords_to_index(n, n - 1, n - 2)]);
    }
}
