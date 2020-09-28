<template>
  <v-container class="text-center">
    <h3>Blog</h3>
    <v-row v-if="allPosts" align="center" justify="center">
      <v-col sm="6">
        <v-expansion-panels class="text-left">
          <v-expansion-panel v-for="(post, idx) in allPosts" :key="idx">
            <v-expansion-panel-header
              expand-icon="mdi-book-open-variant"
              disable-icon-rotate
            >{{ post.title }}</v-expansion-panel-header>
            <v-expansion-panel-content>
              <h5>
                {{
                post.timeOfPost.toDate().toDateString()
                }}
              </h5>
              <br />
              {{ post.post }}
            </v-expansion-panel-content>
          </v-expansion-panel>
        </v-expansion-panels>
      </v-col>
    </v-row>
    <NewPost v-if="auth" />
  </v-container>
</template>

<script>
import { mapGetters, mapActions } from "vuex";
import NewPost from "../components/NewPost.vue";

export default {
  name: "Blog",
  components: {
    NewPost
  },
  methods: {
    ...mapActions(["getPosts", "newLike"]),
    like(post) {
      console.log(post);
      this.newLike(post.id);
    }
  },
  computed: {
    ...mapGetters(["allPosts", "auth"])
  },
  created() {
    this.getPosts().then(console.log("posts fetched"));
  }
};
</script>

<style></style>
