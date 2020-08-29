use piston_window::{Context, rectangle, G2d};
use piston_window::types::Color;

pub const SCALE: f64 = 10.0;

pub fn to_coord(game_coord: i32) -> f64 {
    (game_coord as f64) * SCALE
}

pub fn to_coord_u32(game_coord: i32) -> u32 {
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