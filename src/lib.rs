#![feature(specialization)]

#[macro_use] extern crate pyo3;
#[macro_use] extern crate lazy_static;
extern crate regex;

pub mod alignment;
pub mod record;
pub mod position;

// #[cfg(test)]
// mod tests {
//     #[test]
//     fn it_works() {
//         assert_eq!(2 + 2, 4);
//     }
// }
