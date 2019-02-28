#![feature(specialization)]

#[macro_use] extern crate pyo3;
#[macro_use] extern crate lazy_static;
extern crate regex;

pub mod alignment;
pub mod record;
pub mod position;
pub mod readers;
