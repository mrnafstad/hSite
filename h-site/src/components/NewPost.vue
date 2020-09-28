<template>
  <v-container fluid>
    <v-overlay :value="overlay">
      <v-btn icon @click="overlay = false">
        Post added!
        <v-icon>mdi-check</v-icon>
      </v-btn>
    </v-overlay>
    <v-row justify="center">
      <v-col sm="6">
        <h3>New post to the blog:</h3>
        <v-text-field v-model="post.title" type="string" label="Title..." />
        <v-textarea autoGrow v-model="post.postText" label="Post..." />
        <v-btn class="info" @click="postPost">Submit and publish</v-btn>
      </v-col>
    </v-row>
  </v-container>
</template>

<script>
import { mapActions } from "vuex";
export default {
  name: "NewPost",
  data() {
    return {
      post: {
        postText: "",
        title: ""
      },
      overlay: false
    };
  },
  methods: {
    ...mapActions(["commitPost", "authentication"]),
    postPost() {
      this.commitPost(this.post);
      this.overlay = !this.overlay;
      //this.timer = setTimeout((this.overlay = !this.overlay), 2);
    }
  }
};
</script>

<style>
#newpost {
  border-style: outset;
  background-color: snow;
  padding: 5px;
}

#subPub {
  background-color: lightgrey;
  border-style: outset;
}
</style>
