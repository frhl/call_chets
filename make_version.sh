version="0.1.5"
echo "${version} [$(git rev-parse main | cut -c1-7)]" > .version
