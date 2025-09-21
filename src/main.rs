use ginger::aberth::{aberth, initial_aberth};
use ginger::rootfinding::Options;

fn main() {
    let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
    let mut zrs = initial_aberth(&coeffs);
    let (niter, found) = aberth(&coeffs, &mut zrs, &Options::default());

    if found {
        println!("Found roots in {} iterations:", niter);
        for zr in zrs {
            println!("{}", zr);
        }
    } else {
        println!("Did not find all roots in {} iterations.", niter);
    }
}
