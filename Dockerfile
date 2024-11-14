# Use the latest rocker/verse image (which includes the latest version of R and RStudio)
FROM rocker/verse:latest

# Install the locales package and generate the Icelandic locale
RUN apt-get update && apt-get install -y locales \
    && locale-gen is_IS.UTF-8 \
    && update-locale LANG=is_IS.UTF-8

# Set environment variables for the locale
ENV LANG=is_IS.UTF-8 \
    LANGUAGE=is_IS.UTF-8 \
    LC_ALL=is_IS.UTF-8

# Disable RStudio authentication
ENV DISABLE_AUTH=true

# Expose port 8787 for RStudio
EXPOSE 8787
