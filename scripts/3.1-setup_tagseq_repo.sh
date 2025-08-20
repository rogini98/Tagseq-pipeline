#!/bin/bash

# Utility script to ensure TagSeq repository is properly set up
# Usage: ./setup_tagseq_repo.sh [target_directory]

TARGET_DIR="${1:-/labs/Bolnick/ROL/test-samples}"
REPO_DIR="$TARGET_DIR/tag-based_RNAseq"
REPO_URL="https://github.com/z0on/tag-based_RNAseq.git"

echo "=== TAGSEQ REPOSITORY SETUP ==="
echo "Target directory: $TARGET_DIR"
echo "Repository will be at: $REPO_DIR"
echo ""

# Create target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Check if repository already exists
if [ -d "$REPO_DIR" ]; then
    echo "✅ Repository directory already exists"
    
    # Check if it's actually a git repository
    if [ -d "$REPO_DIR/.git" ]; then
        echo "✅ Valid git repository found"
        
        # Check if tagseq_clipper.pl exists
        if [ -f "$REPO_DIR/tagseq_clipper.pl" ]; then
            echo "✅ tagseq_clipper.pl script found"
            echo ""
            echo "Repository is properly set up and ready to use!"
            exit 0
        else
            echo "❌ tagseq_clipper.pl script missing"
            echo "Repository may be incomplete"
        fi
    else
        echo "⚠️  Directory exists but is not a git repository"
        echo "This may be a partial or corrupted download"
    fi
    
    echo ""
    echo "Options:"
    echo "1. Remove and re-clone: rm -rf $REPO_DIR && $0"
    echo "2. Try updating: cd $REPO_DIR && git pull"
    echo "3. Check manually: ls -la $REPO_DIR"
    
    read -p "Would you like to remove and re-clone? (y/N): " response
    if [[ "$response" =~ ^[Yy]$ ]]; then
        echo "Removing existing directory..."
        rm -rf "$REPO_DIR"
    else
        echo "Keeping existing directory. Please check manually."
        exit 1
    fi
fi

# Clone the repository
echo "Cloning TagSeq repository..."
cd "$TARGET_DIR"

if git clone "$REPO_URL"; then
    echo "✅ Repository cloned successfully"
    
    # Verify the essential files
    if [ -f "$REPO_DIR/tagseq_clipper.pl" ]; then
        echo "✅ tagseq_clipper.pl script found"
        
        # Make the script executable
        chmod +x "$REPO_DIR/tagseq_clipper.pl"
        echo "✅ Made tagseq_clipper.pl executable"
        
        echo ""
        echo "🎉 TagSeq repository setup complete!"
        echo "Repository location: $REPO_DIR"
        echo "Main script: $REPO_DIR/tagseq_clipper.pl"
        
    else
        echo "❌ ERROR: tagseq_clipper.pl not found after cloning"
        echo "Repository may be incomplete or the structure has changed"
        exit 1
    fi
    
else
    echo "❌ ERROR: Failed to clone repository"
    echo ""
    echo "Possible issues:"
    echo "  - No internet connection"
    echo "  - GitHub access blocked"
    echo "  - Repository URL changed"
    echo "  - Insufficient permissions"
    echo ""
    echo "You can try manually:"
    echo "  cd $TARGET_DIR"
    echo "  git clone $REPO_URL"
    exit 1
fi
