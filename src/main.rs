use color_eyre::eyre::Result;
use tracing_subscriber::EnvFilter;

// #![warn(missing_docs)]

#[tokio::main]
async fn main() -> Result<()> {
    setup()?;

    println!("Hello from a (so far completely unnecessary) async runtime");

    Ok(())
}

fn setup() -> Result<()> {
    if std::env::var("RUST_LIB_BACKTRACE").is_err() {
        std::env::set_var("RUST_LIB_BACKTRACE", "1")
    }
    color_eyre::install()?;

    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    tracing_subscriber::fmt::fmt()
        .with_env_filter(EnvFilter::from_default_env())
        .init();

    Ok(())
}
