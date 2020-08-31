use piston_window::types::Color;
use piston_window::*;

mod fluid;
mod draw;

const BACK_COLOR: Color = [0.0, 0.0, 0.0, 1.0];

fn main() {
    let n: usize = 64;

    let mut window: PistonWindow = WindowSettings::new(
        "Fluid simulation",
        [draw::to_coord_u32(n as i32), draw::to_coord_u32(n as i32)],
    )
    .exit_on_esc(true)
    .build()
    .unwrap();

    let mut fluid = fluid::Fluid::new(n, 16, 0.0000001, 0.0001);

    while let Some(event) = window.next() {
        if let Some([mouse_x, mouse_y]) = event.mouse_cursor_args() {
            let fluid_x = (mouse_x / draw::SCALE) as usize;
            let fluid_y = (mouse_y / draw::SCALE) as usize;
            fluid.add_dye(fluid_x, fluid_y, 10.0);
            fluid.add_velocity(fluid_x, fluid_y, 0.0, 1.0)
        }
        window.draw_2d(&event, |c, g, _| {
            clear(BACK_COLOR, g);
            fluid.draw(&c, g);
        });

        event.update(|arg| {
            fluid.update(arg.dt);
        });
    }
}
