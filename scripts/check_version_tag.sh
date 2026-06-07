#!/usr/bin/env bash
# Pre-push guard: if a version tag (v*) is being pushed, assert it
# matches pyproject.toml. Prevents the 1.1.0-on-v1.1.1 class of bug.
#
# Pre-push hooks receive lines on stdin:
#   <local ref> <local sha> <remote ref> <remote sha>

pyproject_version=$(python3 -c "
import tomllib, pathlib
d = tomllib.loads(pathlib.Path('pyproject.toml').read_text())
print(d['project']['version'])
")

while read -r local_ref local_sha remote_ref remote_sha; do
    tag="${local_ref#refs/tags/}"
    case "$tag" in
        v*)
            expected="v${pyproject_version}"
            if [ "$tag" != "$expected" ]; then
                echo "ERROR: pushing tag '$tag' but pyproject.toml says version '${pyproject_version}' (expected tag '$expected')" >&2
                exit 1
            fi
            ;;
    esac
done

exit 0
