# Version Update Instructions

Follow these steps when releasing a new version:

## 1. Update version number
```bash
git tag 1.0.X  # Replace X with your new version number
```

## 2. Generate .version file
```bash
./scripts/make_version.sh
```

## 3. Commit and push
```bash
git add .version
git commit -m "bump version to 1.0.X"
git push origin main
git push origin --tags
```

## Quick one-liner
```bash
git tag 1.0.X && ./scripts/make_version.sh && git add .version && git commit -m "bump version to 1.0.X" && git push origin main && git push origin --tags
```

---

**Note:** Replace `1.0.X` with your actual version number (e.g., `1.0.9`, `1.0.10`, etc.)
