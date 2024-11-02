#![allow(dead_code)]
#![warn(clippy::all)]

// Make sure a graphing calculator agrees with this code.
use std::f64::consts::PI;
use std::cmp::Ordering;
use std::time;
use Func::*;

const NUM_POINTS: usize = 4;
const NUM_COEFFICIENTS: usize = NUM_POINTS - 1;
const BOUNDS: [f64; 2] = [0.0, 0.25 * PI];
const MAX_ADJUSTMENTS: usize = 80;
const ALLOWED_ERROR: f64 = 1E-8;
const EPSILON: f64 = 5E-7;

const REPS: usize = 100;

fn target_function(x: f64) -> f64 {
    x.cos()
}

#[derive(Clone, Debug)]
enum Func {
    Pow(Box<Func>, i32),
    Mul(Box<Func>, Box<Func>),
    Div(Box<Func>, Box<Func>),
    Sub(Box<Func>, Box<Func>),
    Num(f64),
    Var(usize),
}

impl Func {
    fn pack(self) -> Box<Func> {
        Box::new(self)
    }
    fn eval(&self, var: &[f64]) -> f64 {
        match self {
            Func::Pow(a, b) => (a.eval(var)).powi(*b),
            Func::Mul(a, b) => a.eval(var) * b.eval(var),
            Func::Div(a, b) => a.eval(var) / b.eval(var),
            Func::Sub(a, b) => a.eval(var) - b.eval(var),
            Func::Num(n) => *n,
            Func::Var(i) => var[*i],
        }
    }
}

fn main() {
    let mut points = fill_points();
    let start = time::Instant::now();
    let mut total_adjustments = 0;

    println!("{points:?}");

    for _ in 0..REPS {
        points = fill_points();

        for _ in 0..MAX_ADJUSTMENTS {
            total_adjustments += 1;

            let coefficients = coefficients_from_points(&points);
            let extrema = get_extrema(&coefficients, &points);
            let errors: Vec<f64> = extrema
                .iter()
                .map(|extremum| error(&coefficients, *extremum))
                .collect();

            let mut max_absolute_error: f64 = 0.0;
            let mut max_index: usize = 0;

            for (j, error) in errors.iter().enumerate() {
                if (*error).abs() > max_absolute_error {
                    max_absolute_error = (*error).abs();
                    max_index = j;
                }
            }

            let mut min_distance: f64 = f64::MAX;
            let mut min_index: usize = points.len() - 1;

            for (j, point) in points.iter().enumerate() {
                if (error(&coefficients, *point) * errors[max_index]).is_sign_positive()
                    && are_close(*point, extrema[max_index], min_distance)
                {
                    min_distance = (point - extrema[max_index]).abs();
                    min_index = j;
                }
            }

            for i in 0..extrema.len() {
                print!("({}, {}) ", extrema[i], errors[i]);
            }

            println!("\n");

            if are_close(points[min_index], extrema[max_index], ALLOWED_ERROR)
                || (min_index == points.len())
            {
                break;
            }

            points[min_index] = extrema[max_index];
            println!("points {points:?}");
        }
    }

    let time = start.elapsed();
    println!("coefficients {:?}", coefficients_from_points(&points));
    println!("total time {time:?}");
    println!("time per rep {:?}", time / REPS as u32);
    println!("time per adjustment {:?}", time / total_adjustments);
}

// The coefficients of the solved system of equations.
// Go watch a video about the Remez algorithm again.
fn coefficients_from_points(points: &[f64]) -> [f64; NUM_COEFFICIENTS] {
    let mut system = [[0.0; NUM_POINTS + 1]; NUM_POINTS];
    let mut alternating_negative = 1.0;

    // filling values
    for i in 0..NUM_POINTS {
        let mut equation = [0.0; NUM_POINTS + 1];
        let mut point_to_the_n = 1.0;
        alternating_negative = -alternating_negative;

        for coefficient in equation.iter_mut().take(NUM_COEFFICIENTS) {
            *coefficient = point_to_the_n;
            point_to_the_n *= points[i];
            point_to_the_n *= points[i];
        }

        equation[NUM_POINTS - 1] = alternating_negative;
        equation[NUM_POINTS] = target_function(points[i]);

        system[i] = equation;
    }

    for a in 0..NUM_POINTS {
        let product = 1.0 / system[a][a];

        for c in 0..(NUM_POINTS + 1) {
            system[a][c] *= product;
        }

        for b in 0..NUM_POINTS {
            if b == a {
                continue;
            }

            let product = system[b][a];

            for c in 0..(NUM_POINTS + 1) {
                system[b][c] -= system[a][c] * product;
            }
        }
    }

    let mut coefficients: [f64; NUM_COEFFICIENTS] = [0.0; NUM_COEFFICIENTS];

    for i in 0..NUM_COEFFICIENTS {
        coefficients[i] = system[i][NUM_POINTS];
    }

    coefficients
}

fn rref() -> Vec<Func> {
    let mut system: Vec<Vec<Func>> = Vec::with_capacity(NUM_POINTS);
    let mut alternating_negative = 1.0;

    for i in 0..NUM_POINTS {
        let mut equation: Vec<Func> = Vec::with_capacity(NUM_POINTS + 1);
        alternating_negative = -alternating_negative;

        for j in 0..(NUM_POINTS - 1) {
            equation.push(Pow(Var(i).pack(), j as i32));
        }

        equation.push(Num(alternating_negative));
        equation.push(Var(NUM_POINTS + i));

        system.push(equation);
    }

    for a in 0..NUM_POINTS {
        for c in 0..(NUM_POINTS + 1) {
            system[a][c] = match c.cmp(&a) {
                Ordering::Less => Num(0.0),
                Ordering::Equal => Num(1.0),
                Ordering::Greater => Div(system[a][c].clone().pack(), system[a][a].clone().pack())
            }
        }

        for b in 0..NUM_POINTS {
            if b == a { continue; }

            for c in 0..(NUM_POINTS + 1) {
                if a == c {
                    system[b][c] = Num(0.0);
                } else {
                    system[b][c] = Sub(system[b][c].clone().pack(), Mul(system[a][c].clone().pack(),  system[b][a].clone().pack()).pack());
                }
            }
        }
    }

    (0..NUM_COEFFICIENTS).map(|i| system[i][NUM_POINTS].clone()).collect() 
}

fn coefficients_from_points2(points: &[f64], coefficient_functions: &[Func]) -> [f64; NUM_COEFFICIENTS] {
    let mut vars: [f64; NUM_POINTS * 2] = [0.0; NUM_POINTS * 2];

    for i in 0..NUM_POINTS {
        vars[i] = points[i];
        vars[i + NUM_POINTS] = target_function(points[i]);
    }

    let mut coefficients: [f64; NUM_COEFFICIENTS] = [0.0; NUM_COEFFICIENTS];

    for i in 0..NUM_COEFFICIENTS {
        coefficients[i] = coefficient_functions[i].eval(&vars);
    }

    coefficients
}

// The polynomial which will approximate the target function.
fn approximate_polynomial(coefficients: &[f64], x: f64) -> f64 {
    let mut y: f64 = 0.0;
    let mut n = 1.0;

    for coefficient in coefficients {
        y += coefficient * n;
        n *= x;
        n *= x;
    }

    y
}

fn error(coefficients: &[f64], x: f64) -> f64 {
    target_function(x) - approximate_polynomial(coefficients, x)
}

fn get_extrema(coefficients: &[f64], points: &[f64]) -> Vec<f64> {
    let mut extrema: Vec<f64> = Vec::with_capacity(NUM_POINTS * 2);
    extrema.extend_from_slice(&BOUNDS);

    let mut x1: f64 = points[0];
    let mut y1: f64 = error(coefficients, x1);
    let mut slope1: f64 = (error(coefficients, x1 + EPSILON) - y1) / EPSILON;

    for point in points.iter().take(points.len() - 1) {
        // Cubic interpolation using bound points and their slopes.
        let x0: f64 = x1;
        let y0: f64 = y1;
        let slope0: f64 = slope1;

        x1 = *point;
        y1 = error(coefficients, x1);
        slope1 = (error(coefficients, x1 + EPSILON) - y1) / EPSILON;

        // This is impossible to make sense of. Good luck optimizing this without Desmos.
        let width = x0 - x1;
        let height = y0 - y1;
        let width_slope0 = width * slope0;
        let width_slope1 = width * slope1;

        let a = -3.0 * ((-2f64).mul_add(height, width_slope0) + width_slope1);
        let b = (-3f64).mul_add(height, (2f64).mul_add(width_slope0, width_slope1));
        let c = width / a;

        let discriminant = c * (a.mul_add(width_slope0, b * b)).sqrt();
        let the_rest = b.mul_add(c, x0);

        let local_extrema: [f64; 2] = [the_rest + discriminant, the_rest - discriminant];

        for extremum in local_extrema {
            if (BOUNDS[0]..BOUNDS[1]).contains(&extremum) {
                extrema.push(extremum);
            }
        }
    }

    let mut extrema_filtered: Vec<f64> = Vec::with_capacity(extrema.len());

    'outer: for extremum in &extrema {
        for extremum_filtered in &extrema_filtered {
            if are_close(*extremum, *extremum_filtered, ALLOWED_ERROR) {
                continue 'outer;
            }
        }

        extrema_filtered.push(*extremum);
    }

    extrema_filtered
}

fn fill_points() -> [f64; NUM_POINTS] {
    let mut points: [f64; NUM_POINTS] = [0.0; NUM_POINTS];

    for (i, point) in points.iter_mut().enumerate() {
        *point = i as f64 * (BOUNDS[1] - BOUNDS[0]) / (NUM_POINTS - 1) as f64 + BOUNDS[0];
    }

    points
}

fn are_close(num1: f64, num2: f64, error: f64) -> bool {
    (num1 - num2).abs() < error
}
