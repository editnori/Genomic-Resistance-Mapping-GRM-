//! Genomic Resistance Mapping: from BV-BRC records and genome assemblies to
//! interpretable antimicrobial-resistance models and a presentable report.

pub mod analysis;
pub mod annotate;
pub mod bits;
pub mod bvbrc;
pub mod catalog;
pub mod cluster;
pub mod curate;
pub mod dna;
pub mod drugs;
pub mod fasta;
mod fnv;
pub mod kset;
pub mod layout;
pub mod learn;
pub mod matrix;
pub mod panels;
pub mod probes;
pub mod report;
pub mod research;
pub mod serve;
pub mod species;
pub mod species_check;
